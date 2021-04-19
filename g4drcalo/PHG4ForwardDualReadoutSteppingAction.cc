#include "PHG4ForwardDualReadoutSteppingAction.h"
#include "PHG4ForwardDualReadoutDetector.h"

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
PHG4ForwardDualReadoutSteppingAction::PHG4ForwardDualReadoutSteppingAction(PHG4ForwardDualReadoutDetector* detector, const int absorberactive)
  : PHG4SteppingAction(detector->GetName())
  , detector_(detector)
  , hits_(nullptr)
  , absorberhits_(nullptr)
  , hitcontainer(nullptr)
  , hit(nullptr)
  , saveshower(nullptr)
  , absorbertruth(absorberactive)
  , light_scint_model(1)
{
}

PHG4ForwardDualReadoutSteppingAction::~PHG4ForwardDualReadoutSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

//____________________________________________________________________________..
bool PHG4ForwardDualReadoutSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // G4TouchableHistory* theTouchable = (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );  
  G4VPhysicalVolume* volume = touch->GetVolume();
  // G4LogicalVolume* volumeMother = volume->GetLogicalVolume();
  // G4VSensitiveDetector* sensdet = volumeMother->GetSensitiveDetector();
  // G4VSolid* volumeGMother = nullptr;
  // if(volumeMother) volumeGMother = volumeMother->GetSolid();

  // detector_->IsInForwardDualReadout(volume)
  // returns
  //  0 is outside of Forward HCAL
  //  1 is inside scintillator
  // -1 is inside absorber (dead material)

  int whichactive = detector_->IsInForwardDualReadout(volume);


  // if(whichactive){
  // // // if(sensdet)
  // // // cout << theTouchable->GetReplicaNumber() << "\t"<< theTouchable->GetReplicaNumber(1) << "\t"<< theTouchable->GetReplicaNumber(2) << "\t" << volume->GetName() << endl;
  // // // cout << whichactive << "\t" << volume->GetName() << "\tSD:" << sensdet->GetName()  << endl;
  // // // else if(volumeGMother)
  // // // cout << whichactive << "\t" << volume->GetName() << "\tGM:" << volumeGMother->GetName()  << endl;
  // // // else if(volumeMother)
  // // // cout << whichactive << "\t" << volume->GetName() << "\tMM:" << volumeMother->GetName()  << endl;
  // // // else 
  // cout << whichactive << "\t" << volume->GetName() << endl;
  // // // cout << whichactive << "\t" << touch->GetCopyNumber(2)<< "\t" << touch->GetCopyNumber(3) << endl;
  // // // cout << "\tx: " << (aStep->GetPreStepPoint()->GetPosition()).x()<< "\ty: " << (aStep->GetPreStepPoint()->GetPosition()).y()<< "\tz: " << (aStep->GetPreStepPoint()->GetPosition()).z()<< endl;
  // } else {
  // cout << volume->GetCopyNo() << "\t" << volume->GetName() << endl;

  // }


  if (!whichactive)
  {
    return false;
  }

  int layer_id = detector_->get_Layer();
  // unsigned int tower_id = -1;
  int idx_j = -1;
  int idx_k = -1;

  // if (whichactive > 0)  // in sctintillator
  // {
  //   /* Find indizes of sctintillator / tower containing this step */
  //   // FindTowerIndex(touch, idx_j, idx_k);
  //   tower_id = touch->GetVolume(1)->GetCopyNo();
  // }
  // else
  // {
  //   tower_id = touch->GetVolume(1)->GetCopyNo();
  // }

  /* Get energy deposited by this step */
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  G4double light_yield = 0;

  /* Get pointer to associated Geant4 track */
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (detector_->IsBlackHole())
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  /* Make sure we are in a volume */
  if (detector_->IsActive())
  {
    // int idx_l = -1;
    /* Check if particle is 'geantino' */
    bool geantino = false;
    if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
        aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
    {
      geantino = true;
    }

    /* Get Geant4 pre- and post-step points */
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* postPoint = aStep->GetPostStepPoint();
    if(abs( prePoint->GetPosition().y() / cm)>220  || abs( prePoint->GetPosition().x() / cm)>220) return false;
    // cout << "x: " << prePoint->GetPosition().x() << "\ty: "  << prePoint->GetPosition().y() << "\tz: "  << prePoint->GetPosition().z() << endl;
    switch (prePoint->GetStepStatus())
    {
    case fGeomBoundary:
    case fUndefined:
      if (!hit)
      {
        hit = new PHG4Hitv1();
      }
      // hit->set_scint_id(tower_id);
      FindTowerIndexFromPosition(prePoint, idx_j, idx_k);
      // cout << idx_j << ":" << idx_k << endl;
      /* Set hit location (tower index) */
      hit->set_index_j(idx_j);
      hit->set_index_k(idx_k);
      // hit->set_index_l(idx_l);

      // cout << "\tidj: " << idx_j << "\tidk: "  << idx_k << endl;
      // TODO use these positions for an index conversion
      /* Set hit location (space point) */
      hit->set_x(0, prePoint->GetPosition().x() / cm);
      hit->set_y(0, prePoint->GetPosition().y() / cm);
      hit->set_z(0, prePoint->GetPosition().z() / cm);

      /* Set hit time */
      hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);

      //set the track ID
      hit->set_trkid(aTrack->GetTrackID());
      /* set intial energy deposit */
      hit->set_edep(0);
      hit->set_eion(0);
      hit->set_property(PHG4Hit::PROPERTY::scint_gammas, (float)0.);
      hit->set_property(PHG4Hit::PROPERTY::cerenkov_gammas, (float)0.);

      /* Now add the hit to the hit collection */
      if (whichactive > 0)
      {
        hitcontainer = hits_;
        hit->set_light_yield(0);
      }
      else
      {
        hitcontainer = absorberhits_;
      }
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          hit->set_trkid(pp->GetUserTrackId());
          hit->set_shower_id(pp->GetShower()->get_id());
          saveshower = pp->GetShower();
        }
      }
      break;
    default:
      break;
    }

    if (whichactive > 0)
    {
      if (light_scint_model)
      {
        light_yield = GetVisibleEnergyDeposition(aStep);  // for scintillator only, calculate light yields
        static bool once = true;
        if (once && edep > 0)
        {
          once = false;

          if (Verbosity() > 0)
          {
            cout << "PHG4ForwardDualReadoutSteppingAction::UserSteppingAction::"
                 //
                 << detector_->GetName() << " - "
                 << " use scintillating light model at each Geant4 steps. "
                 << "First step: "
                 << "Material = "
                 << aTrack->GetMaterialCutsCouple()->GetMaterial()->GetName()
                 << ", "
                 << "Birk Constant = "
                 << aTrack->GetMaterialCutsCouple()->GetMaterial()->GetIonisation()->GetBirksConstant()
                 << ","
                 << "edep = " << edep << ", "
                 << "eion = " << eion
                 << ", "
                 << "light_yield = " << light_yield << endl;
          }
        }
      }
      else
      {
        light_yield = eion;
      }
    }



    // Double_t fE; //deposited energy in the cell
    // ULong64_t fNphot = 0; // number of optical photons
    float fNscin = 0; // scintillation photons
    float fNcerenkov = 0; // Cerenkov photons

    G4int fScinType; // scintillation process type
    G4int fScinSubType; // scintillation process subtype
    G4int fCerenkovType; // Cerenkov process type
    G4int fCerenkovSubType; // Cerenkov process subtype

    G4Scintillation scin;
    fScinType = scin.GetProcessType();
    fScinSubType = scin.GetProcessSubType();
    G4Cerenkov cer;
    fCerenkovType = cer.GetProcessType();
    fCerenkovSubType = cer.GetProcessSubType();

    //number of optical photons in the event from secondary tracks
    // const std::vector<const G4Track*> *sec = aStep->GetSecondaryInCurrentStep();
    // std::vector<const G4Track*>::const_iterator ittr;
    // for(ittr = sec->begin(); ittr != sec->end(); ittr++) {
      // if((*ittr)->GetParentID() <= 0) continue;

      //all optical photons
    if(aTrack->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
    {
      // fNphot++;

      //identify the process
      G4int ptype = aTrack->GetCreatorProcess()->GetProcessType();
      G4int pstype = aTrack->GetCreatorProcess()->GetProcessSubType();
      // cout << aTrack->GetCreatorProcess()->GetProcessName() << endl;
      //scintillation photons
      G4Material* prevMaterial = aStep->GetPreStepPoint()->GetMaterial();
      if(ptype == fScinType && pstype == fScinSubType && (prevMaterial->GetName().find("G4_POLYSTYRENE") != std::string::npos)){ fNscin++;}
      //Cerenkov photons
      if(ptype == fCerenkovType && pstype == fCerenkovSubType && (prevMaterial->GetName().find("PMMA") != std::string::npos)){ fNcerenkov++;}
    //   if(aTrack->GetParentID() > 0)
    // {
    //   if(aTrack->GetCreatorProcess()->GetProcessName().find("enkov") != string::npos)cout << aTrack->GetCreatorProcess()->GetProcessName() << endl;
    // }
      // if(fNscin>0 || ) cout << aTrack->GetCreatorProcess()->GetProcessName() <<  "\tfNscin: " << fNscin << "\tfNcerenkov: " << fNcerenkov << endl;
      // if(fNscin>0 || fNcerenkov>0) cout << aTrack->GetCreatorProcess()->GetProcessName() <<  "\tfNscin: " << fNscin << "\tfNcerenkov: " << fNcerenkov << endl;
    }//secondary tracks loop
    // cout << __LINE__ << endl;
    //   cout << hit->get_property_float(PHG4Hit::PROPERTY::scint_gammas)<<  "\t" << hit->get_property_float(PHG4Hit::PROPERTY::scint_gammas) << endl;
    // cout << __LINE__ << endl;
    // G4Material* nextMaterial = aStep->GetPostStepPoint()->GetMaterial();
    // string materialstr = nextMaterial->GetName();
    // string materialstr2 = prevMaterial->GetName();
    // if(materialstr.find("G4_AIR") == std::string::npos)cout << materialstr << endl;;
    // if(materialstr2.find("G4_AIR") == std::string::npos)cout << "\t" << materialstr2 << endl;;
      hit->set_property(PHG4Hit::PROPERTY::scint_gammas,hit->get_property_float(PHG4Hit::PROPERTY::scint_gammas)+ fNscin);
    // cout << __LINE__ << endl;
      hit->set_property(PHG4Hit::PROPERTY::cerenkov_gammas,hit->get_property_float(PHG4Hit::PROPERTY::cerenkov_gammas)+ fNcerenkov);
    // cout << __LINE__ << endl;
    //   cout << hit->get_property_float(PHG4Hit::PROPERTY::scint_gammas)<<  "\t" << hit->get_property_float(PHG4Hit::PROPERTY::scint_gammas) << endl;
    // cout << __LINE__ << endl;

    // //number of optical photons in the event from secondary tracks
    // const std::vector<const G4Track*> *sec = aStep->GetSecondaryInCurrentStep();
    // std::vector<const G4Track*>::const_iterator ittr;
    // for(ittr = sec->begin(); ittr != sec->end(); ittr++) {
    //   if((*ittr)->GetParentID() <= 0) continue;

    //   //all optical photons
    //   if((*ittr)->GetDynamicParticle()->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) continue;
    //   fNphot++;

    //   //identify the process
    //   G4int ptype = (*ittr)->GetCreatorProcess()->GetProcessType();
    //   G4int pstype = (*ittr)->GetCreatorProcess()->GetProcessSubType();
    //   //scintillation photons
    //   if(ptype == fScinType && pstype == fScinSubType) fNscin++;
    //   //Cerenkov photons
    //   if(ptype == fCerenkovType && pstype == fCerenkovSubType) fNcerenkov++;

    // }//secondary tracks loop
    // if(sec->size()>0) cout << "fNphot: " << fNphot << "\tfNscin: " << fNscin << "\tfNcerenkov: " << fNcerenkov << endl;
    // }



    /* Update exit values- will be overwritten with every step until
       * we leave the volume or the particle ceases to exist */
    hit->set_x(1, postPoint->GetPosition().x() / cm);
    hit->set_y(1, postPoint->GetPosition().y() / cm);
    hit->set_z(1, postPoint->GetPosition().z() / cm);

    hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    /* sum up the energy to get total deposited */
    hit->set_edep(hit->get_edep() + edep);
    if (whichactive > 0)
    {
      hit->set_eion(hit->get_eion() + eion);
      hit->set_light_yield(hit->get_light_yield() + light_yield);
    }

    if (geantino)
    {
      hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      if (whichactive > 0)
      {
        hit->set_eion(-1);
        hit->set_light_yield(-1);
      }
    }
    if (edep > 0 && (whichactive > 0 || absorbertruth > 0))
    {
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
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
        hitcontainer->AddHit(layer_id, hit);
        if (saveshower)
        {
          saveshower->add_g4hit_id(hitcontainer->GetID(), hit->get_hit_id());
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

    return true;
  }
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4ForwardDualReadoutSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  std::string hitnodename;
  std::string absorbernodename;

  // std::cout << detector_->SuperDetector() << "\t" <<  detector_->GetName() << endl;
  if (detector_->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + detector_->SuperDetector();
    absorbernodename = "G4HIT_ABSORBER_" + detector_->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + detector_->GetName();
    absorbernodename = "G4HIT_ABSORBER_" + detector_->GetName();
  }

  //now look for the map and grab a pointer to it.
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  absorberhits_ = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename);

  // if we do not find the node it's messed up.
  if (!hits_)
  {
    std::cout << "PHG4ForwardDualReadoutSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if (!absorberhits_)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4ForwardDualReadoutSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}

int PHG4ForwardDualReadoutSteppingAction::FindTowerIndexFromPosition(G4StepPoint* prePoint, int& j, int& k)
{
  int j_0 = 0;  //The j and k indices for the scintillator / tower
  int k_0 = 0;  //The j and k indices for the scintillator / tower

  float twrsize = 1.2;//1.2; // was 0.3
  // float twrsize = 1.2;//1.2; // was 0.3
  float drsize = 220.;
  // G4VPhysicalVolume* tower = touch->GetVolume(1);  //Get the tower solid
  // ParseG4VolumeName(tower, j_0, k_0);
  j_0 = (int) ( ( drsize + ( prePoint->GetPosition().x() / cm ) ) / twrsize ); //TODO DRCALO TOWER SIZE
  k_0 = (int) ( ( drsize + ( prePoint->GetPosition().y() / cm ) ) / twrsize ); //TODO DRCALO TOWER SIZE
  // if(prePoint->GetPosition().x()>290) cout << prePoint->GetPosition().x() << "\t" << k_0 << endl;
  j = (j_0 * 1);
  k = (k_0 * 1);

  return 0;
}

int PHG4ForwardDualReadoutSteppingAction::FindTowerIndex(G4TouchableHandle touch, int& j, int& k)
{
  int j_0, k_0;  //The j and k indices for the scintillator / tower

  G4VPhysicalVolume* tower = touch->GetVolume(1);  //Get the tower solid
  ParseG4VolumeName(tower, j_0, k_0);

  j = (j_0 * 1);
  k = (k_0 * 1);

  return 0;
}

int PHG4ForwardDualReadoutSteppingAction::ParseG4VolumeName(G4VPhysicalVolume* volume, int& j, int& k)
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
