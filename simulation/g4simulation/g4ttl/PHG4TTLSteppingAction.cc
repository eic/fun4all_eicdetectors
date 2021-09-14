#include "PHG4TTLSteppingAction.h"
#include "PHG4TTLDetector.h"


class G4VPhysicalVolume;
class PHCompositeNode;
//____________________________________________________________________________..
PHG4TTLSteppingAction::PHG4TTLSteppingAction(PHG4TTLDetector* detector)
  : PHG4SteppingAction(detector->GetName())
  , detector_(detector)
  , hits_(nullptr)
  , hit(nullptr)
  , saveshower(nullptr)
  , layer_id(-1)
  , _isFwd_TTL(false)
  , _N_phi_modules(0)
  , _z_pos_TTL(0)
  , _sensor_resolution_x(0)
  , _sensor_resolution_y(0)
{
}

PHG4TTLSteppingAction::~PHG4TTLSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

//____________________________________________________________________________..
bool PHG4TTLSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  // get volume of the current step
  G4VPhysicalVolume* volume =
      aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;

  const G4Track* aTrack = aStep->GetTrack();

  TVector3 sensorPosition;
  // make sure we are in a volume
  if (detector_->IsInSectorActive(volume))
  {
    bool geantino = false;

    int moduleID = -1;
    int layer = -1;
    int sensor0 = -1;
    int sensor1 = -1;
    int idx_j = -1;
    int idx_k = -1;
    // the check for the pdg code speeds things up, I do not want to make
    // an expensive string compare for every track when we know
    // geantino or chargedgeantino has pid=0
    if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 && aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != std::string::npos)
    {
      geantino = true;
    }
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* postPoint = aStep->GetPostStepPoint();
    //       cout << "track id " << aTrack->GetTrackID() << endl;
    //       cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
    //       cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;
    //layer_id is sector number
    switch (prePoint->GetStepStatus())
    {
    case fGeomBoundary:
    case fUndefined:
      if (!hit)
      {
        hit = new PHG4Hitv1();
      }
      //here we set the entrance values in cm
      hit->set_x(0, prePoint->GetPosition().x() / cm);
      hit->set_y(0, prePoint->GetPosition().y() / cm);
      hit->set_z(0, prePoint->GetPosition().z() / cm);
      // time in ns
      hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
      //set the track ID
      hit->set_trkid(aTrack->GetTrackID());
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          hit->set_trkid(pp->GetUserTrackId());
          hit->set_shower_id(pp->GetShower()->get_id());
          saveshower = pp->GetShower();
        }
      }

      // std::cout << std::endl;
      CalculateSensorHitIndices(prePoint, moduleID,layer,sensor0,sensor1,idx_j, idx_k,sensorPosition);
      hit->set_index_i(moduleID);
      hit->set_index_j(layer);
      hit->set_index_k(sensor0);
      hit->set_index_l(sensor1);
      hit->set_strip_z_index(idx_j);
      hit->set_strip_y_index(idx_k);

      hit->set_local_x(0, sensorPosition.X());
      hit->set_local_y(0, sensorPosition.Y());
      hit->set_local_z(0, sensorPosition.Z());

      //set the initial energy deposit
      hit->set_edep(0);
      hit->set_eion(0);  // only implemented for v5 otherwise empty
      layer_id = 0;//aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(2);
      // std::cout << "layerid: " << layer_id << std::endl;
      //        hit->set_light_yield(0);

      break;
    default:
      break;
    }
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
    hit->set_x(1, postPoint->GetPosition().x() / cm);
    hit->set_y(1, postPoint->GetPosition().y() / cm);
    hit->set_z(1, postPoint->GetPosition().z() / cm);

    hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
    //sum up the energy to get total deposited
    hit->set_edep(hit->get_edep() + edep);
    // std::cout << "energy: " << hit->get_edep() + edep << std::endl;
    hit->set_eion(hit->get_eion() + eion);
    hit->set_path_length(aTrack->GetTrackLength() / cm);
    if (geantino)
    {
      hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
    }
    if (edep > 0)
    {
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp =
                dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
    }
    // if any of these conditions is true this is the last step in
    // this volume and we need to save the hit
    // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
    // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
    // (not sure if this will ever be the case)
    // aTrack->GetTrackStatus() == fStopAndKill: track ends
    if (postPoint->GetStepStatus() == fGeomBoundary ||
        postPoint->GetStepStatus() == fWorldBoundary ||
        postPoint->GetStepStatus() == fAtRestDoItProc ||
        aTrack->GetTrackStatus() == fStopAndKill)
    {
      // save only hits with energy deposit (or -1 for geantino)
      if (hit->get_edep())
      {
        hits_->AddHit(layer_id, hit);
        if (saveshower)
        {
          saveshower->add_g4hit_id(hits_->GetID(), hit->get_hit_id());
        }
        // ownership has been transferred to container, set to null
        // so we will create a new hit for the next track
        hit = nullptr;
      }
      else
      {
        // if this hit has no energy deposit, just reset it for reuse
        // this means we have to delete it in the dtor. If this was
        // the last hit we processed the memory is still allocated
        hit->Reset();
      }
    }

    //       hit->identify();
    // return true to indicate the hit was used
    return true;
  }
  else
  {
    return false;
  }
}


void PHG4TTLSteppingAction::CalculateSensorHitIndices(G4StepPoint* prePoint, int& module_ret, int& layer, int& sensor_0, int& sensor_1, int& j, int& k, TVector3 &sensorposition)
{
  int module_ID = -1;  //The j and k indices for the scintillator / tower
  int layer_ID = -1;  //The j and k indices for the scintillator / tower
  int sensor_ID_0 = -1;  //The j and k indices for the scintillator / tower
  int sensor_ID_1 = -1;  //The j and k indices for the scintillator / tower
  int hit_j_0 = 0;  //The j and k indices for the scintillator / tower
  int hit_k_0 = 0;  //The j and k indices for the scintillator / tower


  G4double baseplate_length = 43.1 * mm;
  G4double baseplate_width = 56.5 * mm / 2;
  G4double baseSH_width = baseplate_width / 2;
  G4double _module_x_dimension = baseplate_length;
  G4double _module_y_dimension = baseplate_width + baseSH_width;

  G4double sensor_y_dimension = 21.2 * mm;
  G4double sensor_x_dimension = 42.0 * mm;

  _sensor_resolution_x = (500e-4 / sqrt(12)) * cm; // in mm
  _sensor_resolution_y = (500e-4 / sqrt(12)) * cm; // in mm

  TVector3 prePointVec(prePoint->GetPosition().x(),prePoint->GetPosition().y(),prePoint->GetPosition().z());
  if(_N_phi_modules>0){
    float prePoint_Phi = prePointVec.Phi()+M_PI;
    module_ID = (int) prePoint_Phi/(2*M_PI/_N_phi_modules);
  }
  if(prePoint->GetPosition().z()>0){
    layer_ID = prePoint->GetPosition().z()>_z_pos_TTL ? 1 : 0;
  } else {
    layer_ID = prePoint->GetPosition().z()<_z_pos_TTL ? 1 : 0;
  }
  if(_isFwd_TTL){
    // front face starts with sensors
    sensor_ID_0 = (int) ( ( ( prePoint->GetPosition().x() + (prePoint->GetPosition().x()<0 ? -_module_x_dimension /2 : _module_x_dimension/2) ) ) / _module_x_dimension );
    sensor_ID_1 = (int) ( ( ( prePoint->GetPosition().y() + (prePoint->GetPosition().y()<0 ? -_module_y_dimension /2 : _module_y_dimension/2) ) ) / _module_y_dimension );
    // calculate x and y position of bottom left corner of sensor
    float sensorcorner_x = sensor_ID_0*_module_x_dimension - sensor_x_dimension/2;
    float sensorcorner_y = 0;
    if((sensor_ID_1%2==0 && layer_ID==0) || (sensor_ID_1%2!=0 && layer_ID==1) ){
    // even sensor counts are located at the bottom of the module
      sensorcorner_y = sensor_ID_1*_module_y_dimension - _module_y_dimension/2 + (0.1 * mm / 2);
    } else {
    // odd sensor counts are located at the top of the module
      sensorcorner_y = sensor_ID_1*_module_y_dimension + _module_y_dimension/2 - sensor_y_dimension - (0.1 * mm / 2);
    }
    // std::cout << "\tcorner_x " << sensorcorner_x << "\tposition_x " <<  prePoint->GetPosition().x() << "\treso " << _sensor_resolution_x << std::endl;
    // std::cout << "\tcorner_y " << sensorcorner_y << "\tposition_y " <<  prePoint->GetPosition().y() << "\treso " << _sensor_resolution_y << std::endl;
    hit_j_0 = (int) ( ( prePoint->GetPosition().x() - sensorcorner_x)  / _sensor_resolution_x );
    sensorposition.SetX( (hit_j_0 * _sensor_resolution_x) + sensorcorner_x );
    hit_k_0 = (int) ( ( prePoint->GetPosition().y() - sensorcorner_y)  / _sensor_resolution_y );
    sensorposition.SetY( (hit_k_0 * _sensor_resolution_y) + sensorcorner_y );

    sensorposition.SetZ( prePoint->GetPosition().z() );
  }


  module_ret = module_ID;
  layer = layer_ID;
  sensor_0 = sensor_ID_0;
  sensor_1 = sensor_ID_1;
  j = hit_j_0;
  k = hit_k_0;
  // if(_isFwd_TTL){
  // std::cout << "module " << module_ID << "\tlayer " << layer_ID << "\tsensor0 " << sensor_ID_0 << "\tsensor1 " << sensor_ID_1 << "\tj " << j << "\tk " << k << std::endl;
  // }
  return;
}


//____________________________________________________________________________..
void PHG4TTLSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  std::string hitnodename;
  if (detector_->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + detector_->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + detector_->GetName();
  }

  //now look for the map and grab a pointer to it.
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());

  // if we do not find the node we need to make it.
  if (!hits_)
  {
    std::cout << "PHG4TTLSteppingAction::SetTopNode - unable to find "
              << hitnodename << std::endl;
  }
}
