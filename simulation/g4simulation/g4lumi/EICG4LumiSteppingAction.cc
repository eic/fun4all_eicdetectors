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

#include "EICG4LumiSteppingAction.h"

#include "EICG4LumiDetector.h"
#include "EICG4LumiSubsystem.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4StepStatusDecode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

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
EICG4LumiSteppingAction::EICG4LumiSteppingAction(EICG4LumiSubsystem *subsys, EICG4LumiDetector *detector, const PHParameters *parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Subsystem(subsys)
  , m_Detector(detector)
  , m_Params(parameters)
  , m_HitContainerCAL(nullptr)
  , m_HitContainerTracking(nullptr)
  , m_HitContainerVirt(nullptr)
  , m_Hit(nullptr)
  , m_SaveShower(nullptr)
  , m_SaveVolPre(nullptr)
  , m_SaveVolPost(nullptr)
  , m_SaveLightYieldFlag(m_Params->get_int_param("lightyield"))
  , m_SaveTrackId(-1)
  , m_SavePreStepStatus(-1)
  , m_SavePostStepStatus(-1)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_BlackHoleFlag(m_Params->get_int_param("blackhole"))
  , m_UseG4StepsFlag(m_Params->get_int_param("use_g4steps"))
  , m_Tmin(m_Params->get_double_param("tmin") * ns)
  , m_Tmax(m_Params->get_double_param("tmax") * ns)
  , m_EdepSum(0)
  , m_EabsSum(0)
  , m_EionSum(0)
{
// G4 seems to have issues in the um range
   m_Zmin -= copysign(m_Zmin, 1. / 1e6 * cm);
   m_Zmax += copysign(m_Zmax, 1. / 1e6 * cm);
}

//____________________________________________________________________________..
EICG4LumiSteppingAction::~EICG4LumiSteppingAction()
{
   delete m_Hit;
}

//____________________________________________________________________________..
// This is the implementation of the G4 UserSteppingAction
bool EICG4LumiSteppingAction::UserSteppingAction(const G4Step *aStep,bool was_used)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  G4double light_yield = GetVisibleEnergyDeposition(aStep);

  // get volume of the current step
  G4VPhysicalVolume *volume = touch->GetVolume();
  G4LogicalVolume *volumeLogical = volume->GetLogicalVolume();
  G4Material *mat = volumeLogical->GetMaterial();

  int activeMaterial = m_Detector->IsInDetector(volume);
  int virtualMaterial = m_Detector->IsInVirtualDetector(volume);
  int trackingMaterial = (mat->GetName().compareTo("G4_Si")==0) ? 1 : 0;

  if ( !activeMaterial && !virtualMaterial )
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

  int layer_id = m_Detector->get_Layer();

  int layer_type = 1;
  if (!m_ActiveFlag)
  {
	return false;
  }
  bool geantino = false;

  if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
      aTrack->GetParticleDefinition()->GetParticleName().find("geantino") !=
          std::string::npos)  // this also accounts for "chargedgeantino"
  {
    geantino = true;
  }
  G4StepPoint *prePoint = aStep->GetPreStepPoint();
  G4StepPoint *postPoint = aStep->GetPostStepPoint();

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
// These are the normal cases

  case fGeomBoundary:
  case fUndefined:
    if (!m_Hit)
    {
      m_Hit = new PHG4Hitv1();
    }
    m_Hit->set_layer((unsigned int) layer_id);
    // here we set the entrance values in cm
    //FindTowerIndexFromPosition(prePoint, idx_j, idx_k);
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

    if (m_SaveLightYieldFlag)
    {
       		m_Hit->set_light_yield(0);
    }

   if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
      {
        m_Hit->set_trkid(pp->GetUserTrackId());
	m_Hit->set_shower_id(pp->GetShower()->get_id());
	m_SaveShower = pp->GetShower();
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
    m_Hit->set_eion(-1);
    m_Hit->set_light_yield(-1);
	
  }
  else
  {
    if (!m_BlackHoleFlag)
    {
       m_Hit->set_eion(m_Hit->get_eion() + eion);
  if(m_SaveLightYieldFlag)
	  	m_Hit->set_light_yield(m_Hit->get_light_yield()+light_yield);
    }
  }

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
      if( activeMaterial ) 
      {
        if( trackingMaterial ) {
          m_HitContainerTracking->AddHit(layer_id, m_Hit);
        }
        else {
          m_HitContainerCAL->AddHit(layer_id, m_Hit);
        }
      }
      if( virtualMaterial ) 
      {
	m_HitContainerVirt->AddHit(layer_id, m_Hit);
      }

      if (m_SaveShower)
      {
          m_SaveShower->add_g4hit_id(m_HitContainerCAL->GetID(), m_Hit->get_hit_id());
      }
      m_Hit = nullptr;
    }
    else
    {
           m_Hit->Reset();
    }
  }

 // return true to indicate the hit was used
  return true;
}

//____________________________________________________________________________..
void EICG4LumiSteppingAction::SetInterfacePointers(PHCompositeNode *topNode)
{
  //std::string hitnodename = "G4HIT_" + m_Detector->GetName();

  // now look for the map and grab a pointer to it.
  m_HitContainerCAL = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeNameCAL);
  m_HitContainerTracking = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeNameTracking);
  m_HitContainerVirt = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeNameVirt);

  // if we do not find the node we need to make it.
  if (!m_HitContainerCAL)
  {
    std::cout << PHWHERE << " EICG4LumiSteppingAction::SetTopNode - unable to find "
              << m_HitNodeNameCAL << std::endl;
    
    gSystem->Exit(1);
  }
  if (!m_HitContainerTracking)
  {
    std::cout << PHWHERE << " EICG4LumiSteppingAction::SetTopNode - unable to find "
              << m_HitNodeNameTracking << std::endl;
    
    gSystem->Exit(1);
  }
  if (!m_HitContainerVirt)
  {
    std::cout << PHWHERE << " EICG4LumiSteppingAction::SetTopNode - unable to find "
              << m_HitNodeNameVirt << std::endl;
    
    gSystem->Exit(1);
  }
}

//__________________________________________________________________________..
bool EICG4LumiSteppingAction::hasMotherSubsystem() const
{
  if (m_Subsystem->GetMotherSubsystem())
  {
    return true;
  }
  return false;
}
