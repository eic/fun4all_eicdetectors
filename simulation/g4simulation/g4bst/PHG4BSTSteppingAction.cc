#include "PHG4BSTSteppingAction.h"
#include "PHG4BSTDetector.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>         // for PHG4SteppingAction

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4IonisParamMat.hh>           // for G4IonisParamMat
#include <Geant4/G4Material.hh>                // for G4Material
#include <Geant4/G4MaterialCutsCouple.hh>
#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepStatus.hh>              // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>                  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>             // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>         // for G4TouchableHandle
#include <Geant4/G4Track.hh>                   // for G4Track
#include <Geant4/G4VSolid.hh>                   // for G4Track
#include <Geant4/G4TrackStatus.hh>             // for fStopAndKill
#include <Geant4/G4Types.hh>                   // for G4double
#include <Geant4/G4VPhysicalVolume.hh>         // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>              // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>   // for G4VUserTrackInformation
#include <Geant4/G4VSensitiveDetector.hh>   // for G4VUserTrackInformation
#include <Geant4/G4OpticalPhoton.hh>
#include <Geant4/G4Scintillation.hh>
#include <Geant4/G4Cerenkov.hh>
#include <TSystem.h>

#include <boost/tokenizer.hpp>
// this is an ugly hack, the gcc optimizer has a bug which
// triggers the uninitialized variable warning which
// stops compilation because of our -Werror
#include <boost/version.hpp>  // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 4 && BOOST_VERSION == 105700)
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma message "ignoring bogus gcc warning in boost header lexical_cast.hpp"
#include <boost/lexical_cast.hpp>
#pragma GCC diagnostic warning "-Wuninitialized"
#else
#include <boost/lexical_cast.hpp>
#endif

#include <iostream>
#include <string>                              // for basic_string, operator+

class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
PHG4BSTSteppingAction::PHG4BSTSteppingAction(PHG4BSTDetector* detector, const int absorberactive)
  : PHG4SteppingAction(detector->GetName())
  , detector_(detector)
  , hits_(0)
  , absorberhits_(nullptr)
  , hitcontainer(nullptr)
  , hit(nullptr)
  , saveshower(nullptr)
  , _towerdivision(0.0)
  , _tower_size(1.0)
  , _readout_size(1.0)
  , _detector_size(100)
  , absorbertruth(absorberactive)
  , light_scint_model(1)
{
  hits_ = new PHG4HitContainer*[6];
}

PHG4BSTSteppingAction::~PHG4BSTSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

//____________________________________________________________________________..
bool PHG4BSTSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  // get volume of the current step
  G4VPhysicalVolume* volume =
      aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;

  const G4Track* aTrack = aStep->GetTrack();


  if (detector_->IsInActiveSensorBST(volume))
  {
    
    int layer_id = -1;
    std::string bstLayerNameFind = "BST_";
    if (volume->GetName().find(bstLayerNameFind) != std::string::npos) {
      auto pos = volume->GetName().find(bstLayerNameFind);
      layer_id = std::stoi(volume->GetName().substr(pos + bstLayerNameFind.size(), pos + bstLayerNameFind.size() + 1));
    }
    // cout << volume->GetName() << "\t" << layer_id << "\t" << layer_id << endl;
    if(layer_id<0 || layer_id>6){
      cout << "ERROR: BST layer id is out of range" << endl;
      return false;
    }

    bool geantino = false;

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
      hit->set_index_i(layer_id);

      //set the initial energy deposit
      hit->set_edep(0);
      hit->set_eion(0);  // only implemented for v5 otherwise empty
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
        hits_[layer_id]->AddHit(layer_id, hit);
        if (saveshower)
        {
          saveshower->add_g4hit_id(hits_[layer_id]->GetID(), hit->get_hit_id());
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

//____________________________________________________________________________..
void PHG4BSTSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  std::string hitnodename;
  std::string absorbernodename;

  // std::cout << detector_->SuperDetector() << "\t" <<  detector_->GetName() << endl;
  if (detector_->SuperDetector() != "NONE")
  {
    absorbernodename = "G4HIT_ABSORBER_" + detector_->SuperDetector();
  }
  else
  {
    absorbernodename = "G4HIT_ABSORBER_" + detector_->GetName();
  }

  for(int ilay=0; ilay<6; ilay++){
    if (detector_->SuperDetector() != "NONE")
    {
      hitnodename = "G4HIT_" + detector_->SuperDetector() + "_" + std::to_string(ilay);
    }
    else
    {
      hitnodename = "G4HIT_" + detector_->GetName() + "_" + std::to_string(ilay);
    }
    //now look for the map and grab a pointer to it.
    hits_[ilay] = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);

    // if we do not find the node it's messed up.
    if (!hits_[ilay])
    {
      std::cout << "PHG4BSTSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
      gSystem->Exit(1);
  }
  }
  absorberhits_ = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename);
  if (!absorberhits_)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4BSTSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}

int PHG4BSTSteppingAction::ParseG4VolumeName(G4VPhysicalVolume* volume, int& j, int& k)
{
  // cout << volume->GetName() << endl;
  boost::char_separator<char> sep("_");
  boost::tokenizer<boost::char_separator<char> > tok(volume->GetName(), sep);
  boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;
  for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
  {
    if (*tokeniter == "j")
    {
      ++tokeniter;
      if (tokeniter == tok.end()) break;
      j = boost::lexical_cast<int>(*tokeniter);
    }
    else if (*tokeniter == "k")
    {
      ++tokeniter;
      if (tokeniter == tok.end()) break;
      k = boost::lexical_cast<int>(*tokeniter);
    }
  }
  return 0;
}
