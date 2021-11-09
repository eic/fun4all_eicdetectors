//____________________________________________________________________________..
//
// This is a working template for the Stepping Action which needs to be implemented
// for active detectors. Most of the code is error handling and access to the G4 objects
// and our data structures. It does not need any adjustment. The only thing you need to
// do is to add the properties of the G4Hits you want to save for later analysis
// This needs to be done in 2 places, G4Hits are generated when a G4 track enters a new
// volume (or is created). Here you give it an initial value. When the G4 track leaves
// the volume the final value needs to be set.
// The places to do this is marked by //implement your own here//
//
// As guidance you can look at the total (integrated over all steps in a volume) energy
// deposit which should always be saved.
// Additionally the total ionization energy is saved - this can be removed if you are not
// interested in this. Naturally you may want remove these comments in your version
//
//____________________________________________________________________________..

#include "EICG4B0SteppingAction.h"

#include "EICG4B0Detector.h"
#include "EICG4B0Subsystem.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4StepStatusDecode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4NavigationHistory.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ReferenceCountedHandle.hh>
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>
#include <Geant4/G4StepStatus.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4TrackStatus.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VTouchable.hh>
#include <Geant4/G4VUserTrackInformation.hh>

#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>

class PHCompositeNode;

//____________________________________________________________________________..
EICG4B0SteppingAction::EICG4B0SteppingAction(EICG4B0Subsystem *subsys, EICG4B0Detector *detector, const PHParameters *parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Subsystem(subsys)
  , m_Detector(detector)
  , m_Params(parameters)
  , m_HitContainer(nullptr)
  , m_Hit(nullptr)
  , m_SaveShower(nullptr)
  , m_SaveVolPre(nullptr)
  , m_SaveVolPost(nullptr)
//  , m_SaveLightYieldFlag(m_Params->get_int_param("lightyield"))
  , m_SaveTrackId(-1)
  , m_SavePreStepStatus(-1)
  , m_SavePostStepStatus(-1)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_BlackHoleFlag(m_Params->get_int_param("blackhole"))
  , m_UseG4StepsFlag(m_Params->get_int_param("use_g4steps"))
  , m_Zmin(m_Params->get_double_param("place_z") * cm - m_Params->get_double_param("length") * cm / 2.)
  , m_Zmax(m_Params->get_double_param("place_z") * cm + m_Params->get_double_param("length") * cm / 2.)
  , m_Tmin(m_Params->get_double_param("tmin") * ns)
  , m_Tmax(m_Params->get_double_param("tmax") * ns)
  , m_EdepSum(0)
  , m_EabsSum(0)
  , _towerdivision(0.0)
  , _tower_size(2.0)
  , _readout_size(2.0)
  , _detector_size(20)
{
// G4 seems to have issues in the um range
   m_Zmin -= copysign(m_Zmin, 1. / 1e6 * cm);
   m_Zmax += copysign(m_Zmax, 1. / 1e6 * cm);
}

//____________________________________________________________________________..
EICG4B0SteppingAction::~EICG4B0SteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
// This is the implementation of the G4 UserSteppingAction
bool EICG4B0SteppingAction::UserSteppingAction(const G4Step *aStep, bool was_used)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume *volume = touch->GetVolume();
  // IsInDetector(volume) returns
  //  == 0 outside of detector
  //   > 0 for hits in active volume
  //  < 0 for hits in passive material
  int whichactive = m_Detector->IsInDetector(volume);
  if (!whichactive)
  {
    return false;
  }
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  const G4Track *aTrack = aStep->GetTrack();
  // if this detector stops everything, just put all kinetic energy into edep
  if (m_BlackHoleFlag)
  {
    if ((!std::isfinite(m_Tmin) && !std::isfinite(m_Tmax)) ||
          aTrack->GetGlobalTime() < m_Tmin ||
          aTrack->GetGlobalTime() > m_Tmax)
	{
	    edep = aTrack->GetKineticEnergy() / GeV;
	    G4Track *killtrack = const_cast<G4Track *>(aTrack);
	    killtrack->SetTrackStatus(fStopAndKill);
	}
  }
  // we use here only one detector in this simple example
  // if you deal with multiple detectors in this stepping action
  // the detector id can be used to distinguish between them
  // hits can easily be analyzed later according to their detector id
  int layer_id = m_Detector->get_Layer();
  int idx_j = -1;
  int idx_k = -1;
// deadhits for dead material
 int layer_type = 0;
  if (m_Params->get_string_param("material") == "G4_Cu") layer_type=0;
  else layer_type = 1; 
  //if (layer_id % 2 == 1) //Proper implementation for several layer configurations
/*
   if (m_Params->get_string_param("material") == "G4_Cu") 
  {
	return false;
  }
*/
  if (!m_ActiveFlag)
  {
	return false;
  }
  bool geantino = false;
  // the check for the pdg code speeds things up, I do not want to make
  // an expensive string compare for every track when we know
  // geantino or chargedgeantino has pid=0
  if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
      aTrack->GetParticleDefinition()->GetParticleName().find("geantino") !=
          std::string::npos)  // this also accounts for "chargedgeantino"
  {
    geantino = true;
  }
  G4StepPoint *prePoint = aStep->GetPreStepPoint();
  G4StepPoint *postPoint = aStep->GetPostStepPoint();

  // Here we have to decide if we need to create a new hit.  Normally this should
  // only be neccessary if a G4 Track enters a new volume or is freshly created
  // For this we look at the step status of the prePoint (beginning of the G4 Step).
  // This should be either fGeomBoundary (G4 Track crosses into volume) or
  // fUndefined (G4 Track newly created)
  // Sadly over the years with different G4 versions we have observed cases where
  // G4 produces "impossible hits" which we try to catch here
  // These errors were always rare and it is not clear if they still exist but we
  // still check for them for safety. We can reproduce G4 runs identically (if given
  // the sequence of random number seeds you find in the log), the printouts help
  // us giving the G4 support information about those failures
  //
  switch (prePoint->GetStepStatus())
  {
  case fPostStepDoItProc:
    if (m_SavePostStepStatus != fGeomBoundary)
    {
      // this is the okay case, fPostStepDoItProc called in a volume, not first thing inside
      // a new volume, just proceed here
      break;
    }
    else
    {
      // this is an impossible G4 Step print out diagnostic to help debug, not sure if
      // this is still with us
      std::cout << GetName() << ": New Hit for  " << std::endl;
      std::cout << "prestep status: "
                << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
                << ", poststep status: "
                << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
                << ", last pre step status: "
                << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
                << ", last post step status: "
                << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
      std::cout << "last track: " << m_SaveTrackId
                << ", current trackid: " << aTrack->GetTrackID() << std::endl;
      std::cout << "phys pre vol: " << volume->GetName()
                << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
      std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
                << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
    }
    break;
    // These are the normal cases
  case fGeomBoundary:
  case fUndefined:
    if (!m_Hit)
    {
      m_Hit = new PHG4Hitv1();
    }
    m_Hit->set_layer((unsigned int) layer_id);
    // here we set the entrance values in cm
    FindTowerIndexFromPosition(prePoint, idx_j, idx_k);
  //std::cout << "B0 Hits: "<<postPoint->GetPosition().x() / cm<<" "<<postPoint->GetPosition().y() / cm<<" "<<edep<<std::endl;

    m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
    m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
    m_Hit->set_z(0, prePoint->GetPosition().z() / cm);
   
    m_Hit->set_px(0, prePoint->GetMomentum().x() / GeV);
    m_Hit->set_py(0, prePoint->GetMomentum().y() / GeV);
    m_Hit->set_pz(0, prePoint->GetMomentum().z() / GeV);
  
    // time in ns
    m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
    // set the track ID
    m_Hit->set_trkid(aTrack->GetTrackID());
    m_SaveTrackId = aTrack->GetTrackID();
    // set the initial energy deposit
    m_EdepSum = 0;
    m_EabsSum = 0;

    m_Hit->set_edep(0);
    if (!geantino && !m_BlackHoleFlag)
    {
           m_Hit->set_eion(0);
    }
    if (layer_type){
//	m_SaveHitContainer=m_HitContainer;
	m_Hit->set_eion(0);
//	if (m_SaveLightYieldFlag)
//  	{
//       		m_Hit->set_light_yield(0);
		m_Hit->set_index_j(idx_j);
	    	m_Hit->set_index_k(idx_k);
//   	}
    }
/*    else{
	m_SaveHitContainer=m_AbsorberHitContainer;
   }*/
    // implement your own here://
    // add the properties you are interested in via set_XXX methods
    // you can find existing set methods in $OFFLINE_MAIN/include/g4main/PHG4Hit.h
    // this is initialization of your value. This is not needed you can just set the final
    // value at the last step in this volume later one
/*    if (whichactive > 0)
    {
      m_EionSum = 0;  // assuming the ionization energy is only needed for active
                      // volumes (scintillators)
      m_Hit->set_eion(0);
      m_SaveHitContainer = m_HitContainer;
    }
    else
    {
      std::cout << "implement stuff for whichactive < 0 (inactive volumes)" << std::endl;
      gSystem->Exit(1);
    }
*/    // this is for the tracking of the truth info
    if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
      {
        m_Hit->set_trkid(pp->GetUserTrackId());
	m_Hit->set_shower_id(pp->GetShower()->get_id());
	m_SaveShower = pp->GetShower();
//        pp->GetShower()->add_g4hit_id(m_SaveHitContainer->GetID(), m_Hit->get_hit_id());
      }
    }
    if (!hasMotherSubsystem() && (m_Hit->get_z(0) * cm > m_Zmax || m_Hit->get_z(0) * cm < m_Zmin))
    {
        std::cout << m_Detector->SuperDetector() << std::setprecision(9)
        << " PHG4CylinderSteppingAction: Entry hit z " << m_Hit->get_z(0) * cm
        << " outside acceptance,  zmin " << m_Zmin
        << ", zmax " << m_Zmax << ", layer: " << layer_id << std::endl;
    }
    break;
  default:
    break;
  }

  // This section is called for every step
  // some sanity checks for inconsistencies (aka bugs) we have seen over the years
  // check if this hit was created, if not print out last post step status
  if (!m_Hit || !std::isfinite(m_Hit->get_x(0)))
  {
    std::cout << GetName() << ": hit was not created" << std::endl;
    std::cout << "prestep status: "
              << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
              << ", poststep status: "
              << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
              << ", last pre step status: "
              << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
              << ", last post step status: "
              << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
    std::cout << "last track: " << m_SaveTrackId
              << ", current trackid: " << aTrack->GetTrackID() << std::endl;
    std::cout << "phys pre vol: " << volume->GetName()
              << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
    std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
              << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
    // This is fatal - a hit from nowhere. This needs to be looked at and fixed
    gSystem->Exit(1);
  }
  m_SavePostStepStatus = postPoint->GetStepStatus();
  // check if track id matches the initial one when the hit was created
  if (aTrack->GetTrackID() != m_SaveTrackId)
  {
    std::cout << GetName() << ": hits do not belong to the same track" << std::endl;
    std::cout << "saved track: " << m_SaveTrackId
              << ", current trackid: " << aTrack->GetTrackID()
              << ", prestep status: " << prePoint->GetStepStatus()
              << ", previous post step status: " << m_SavePostStepStatus << std::endl;
    // This is fatal - a hit from nowhere. This needs to be looked at and fixed
    gSystem->Exit(1);
  }

  // We need to cache a few things from one step to the next
  // to identify impossible hits and subsequent debugging printout
  m_SavePreStepStatus = prePoint->GetStepStatus();
  m_SavePostStepStatus = postPoint->GetStepStatus();
  m_SaveVolPre = volume;
  m_SaveVolPost = touchpost->GetVolume();
  
  m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
  m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
  m_Hit->set_z(1, postPoint->GetPosition().z() / cm);
  
  m_Hit->set_px(1, postPoint->GetMomentum().x() / GeV);
  m_Hit->set_py(1, postPoint->GetMomentum().y() / GeV);
  m_Hit->set_pz(1, postPoint->GetMomentum().z() / GeV);
  
  m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
  //sum up the energy to get total deposited
  m_Hit->set_edep(m_Hit->get_edep() + edep);
  if (layer_type){
	m_Hit->set_eion(m_Hit->get_eion()+eion);
	m_Hit->set_light_yield(m_Hit->get_light_yield()+eion);
  }
  if (!layer_type) m_EabsSum += edep; 
  m_EdepSum += edep;
  if (!hasMotherSubsystem() && (m_Hit->get_z(1) * cm > m_Zmax || m_Hit->get_z(1) * cm < m_Zmin))
  {
    std::cout << m_Detector->SuperDetector() << std::setprecision(9)
        << " PHG4CylinderSteppingAction: Exit hit z " << m_Hit->get_z(1) * cm
        << " outside acceptance zmin " << m_Zmin
        << ", zmax " << m_Zmax << ", layer: " << layer_id << std::endl;
  }
  if (geantino)
  {
    m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
	if(layer_type){
		m_Hit->set_eion(-1);
//		m_Hit->set_light_yield(-1);
	}
  }
  else
  {
    if (!m_BlackHoleFlag)
    {
       double eion = edep - aStep->GetNonIonizingEnergyDeposit() / GeV;
       m_Hit->set_eion(m_Hit->get_eion() + eion);
    }
  }

/*  if (m_SaveLightYieldFlag)
  {
    double light_yield = GetVisibleEnergyDeposition(aStep) / GeV;
    m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield);
  }
*/
  if (edep > 0 || m_SaveAllHitsFlag)
  {
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
       {
         pp->SetKeep(1);  // we want to keep the track
       }
    }
  }
  // here we just update the exit values, it will be overwritten
  // for every step until we leave the volume or the particle
  // ceases to exist
  // sum up the energy to get total deposited
  ///m_EdepSum += edep;
  ///if (whichactive > 0)
  ///{
  ///  m_EionSum += eion;
  ///}
  // if any of these conditions is true this is the last step in
  // this volume and we need to save the hit
  // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
  // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
  // (happens when your detector goes outside world volume)
  // postPoint->GetStepStatus() == fAtRestDoItProc: track stops (typically
  // aTrack->GetTrackStatus() == fStopAndKill is also set)
  // aTrack->GetTrackStatus() == fStopAndKill: track ends
  if (postPoint->GetStepStatus() == fGeomBoundary ||
      postPoint->GetStepStatus() == fWorldBoundary ||
      postPoint->GetStepStatus() == fAtRestDoItProc ||
      aTrack->GetTrackStatus() == fStopAndKill ||
      m_UseG4StepsFlag > 0)
  {
    // save only hits with energy deposit (or geantino)
    if (m_Hit->get_edep() || m_SaveAllHitsFlag)
    {
	m_Hit->set_layer(layer_id);
	m_Hit->set_hit_type(layer_type);
//	m_Hit->set_eion(m_EionSum);
//	m_Hit->set_edep(m_EdepSum);
      m_HitContainer->AddHit(layer_id, m_Hit);
      if (m_SaveShower)
      {
          m_SaveShower->add_g4hit_id(m_HitContainer->GetID(), m_Hit->get_hit_id());
      }
      m_Hit = nullptr;
    }
    else
    {
      // if this hit has no energy deposit, just reset it for reuse
      // this means we have to delete it in the dtor. If this was
      // the last hit we processed the memory is still allocated
      m_Hit->Reset();
    }
  }
  // return true to indicate the hit was used
  return true;
}

//____________________________________________________________________________..
int EICG4B0SteppingAction::FindTowerIndexFromPosition(G4StepPoint* prePoint, int& j, int& k)
{
  int j_0 = 0;  //The j and k indices for the scintillator / tower
  int k_0 = 0;  //The j and k indices for the scintillator / tower
  int towersize = 2;
  if(_towerdivision==0.){
    int maxsubtow = (int) ( (_tower_size) / (_readout_size));
    _towerdivision = (_tower_size - (maxsubtow * _readout_size))/maxsubtow;
    _towerdivision+=_readout_size;
  }
//  j_0 = (int) ( ( _detector_size + ( prePoint->GetPosition().x() ) ) / _towerdivision ); //TODO DRCALO TOWER SIZE
//  k_0 = (int) ( ( _detector_size + ( prePoint->GetPosition().y() ) ) / _towerdivision ); //TODO DRCALO TOWER SIZE

  j_0 = (int) ( ( 36.6 + ( prePoint->GetPosition().x() / cm ) ) / towersize ); //TODO DRCALO TOWER SIZE
  k_0 = (int) ( ( 22 + ( prePoint->GetPosition().y() / cm ) ) / towersize ); //TODO DRCALO TOWER SIZE
  j = (j_0 * 1);
  k = (k_0 * 1);

  return 0;
}
//____________________________________________________________________________..
void EICG4B0SteppingAction::SetInterfacePointers(PHCompositeNode *topNode)
{
  //std::string hitnodename = "G4HIT_" + m_Detector->GetName();
  //  std::cout << " ---> !!! hitnodename: " << hitnodename << std::endl;
  // now look for the map and grab a pointer to it.
//  std::string hitnodename;
//  std::string absorbernodename;
//  if (m_Detector->SuperDetector()!= "NONE")
//  {
//	hitnodename = "G4HIT_" + m_Detector->SuperDetector();
//	absorbernodename = "G4HIT_ABSORBER_" + m_Detector->SuperDetector();
//  }
//  else
//  {
//	 hitnodename = "G4HIT_" + m_Detector->GetName();
//	 absorbernodename = "G4HIT_ABSORBER_" + m_Detector->GetName();
//  }
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
//  m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename);
  // if we do not find the node we need to make it.
  if (!m_HitContainer)
  {
    std::cout << "EICG4B0SteppingAction::SetTopNode - unable to find "
              << m_HitNodeName << std::endl;
	gSystem->Exit(1);
  }
/* if (!m_AbsorberHitContainer)
  {
    if (Verbosity() > 0)
    {
      std::cout << "EICG4B0SteppingAction::SetInterfacePointers - unable to find " << absorbernodename << std::endl;
	std::cout <<"Are you running Realistic B0 ???"<<std::endl;
    }
  }*/
}
bool EICG4B0SteppingAction::hasMotherSubsystem() const
{
  if (m_Subsystem->GetMotherSubsystem())
  {
    return true;
  }
  return false;
}
