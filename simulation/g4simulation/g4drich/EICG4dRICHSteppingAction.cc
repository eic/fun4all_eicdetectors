#include "EICG4dRICHSteppingAction.h"
#include "EICG4dRICHDetector.h"
#include "EICG4dRICHHit.h"

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllBase.h>

#include <g4detectors/PHG4StepStatusDecode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <TSystem.h>

#include <G4ParticleDefinition.hh>
#include <G4ReferenceCountedHandle.hh>
#include <G4Step.hh>
#include <G4StepStatus.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4TouchableHandle.hh>
#include <G4TrackStatus.hh>
#include <G4Types.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VProcess.hh>
#include <G4VTouchable.hh>
#include <G4VUserTrackInformation.hh>

#include <cmath>
#include <iostream>
#include <string>

class PHCompositeNode;

//____________________________________________________________________________..
EICG4dRICHSteppingAction::EICG4dRICHSteppingAction(EICG4dRICHDetector *detector,
                                         const PHParameters *parameters)
    : PHG4SteppingAction(detector->GetName())
    , m_Detector(detector)
    , m_Params(parameters)
    , m_HitContainer(nullptr)
    , m_Hit(nullptr)
    , m_SaveHitContainer(nullptr)
    , m_SaveVolPre(nullptr)
    , m_SaveVolPost(nullptr)
    , m_SaveTrackId(-1)
    , m_SavePreStepStatus(-1)
    , m_SavePostStepStatus(-1)
    , m_ActiveFlag(m_Params->get_int_param("active"))
    , m_EdepSum(0)
    , m_EionSum(0)
{
  // hit type strings
  hitTypeStr[hEntrance] = "entrance";
  hitTypeStr[hExit] = "exit";
  hitTypeStr[hPSST] = "psst";
  hitTypeStr[hIgnore] = "ignore";
  // hit subtype strings
  // - entrances
  hitSubtypeStr[entPrimary] = "primary";
  hitSubtypeStr[entSecondary] = "secondary";
  hitSubtypeStr[entPostStep] = "postStep";
  // - exits
  hitSubtypeStr[exPrimary] = "primary";
  hitSubtypeStr[exSecondary] = "secondary";
  // - photosensor hits
  hitSubtypeStr[psOptical] = "optical";
  hitSubtypeStr[psGamma] = "gamma";
  hitSubtypeStr[psOther] = "other";
  // - unknown
  hitSubtypeStr[subtypeUnknown] = "unknown";
}

//____________________________________________________________________________..
EICG4dRICHSteppingAction::~EICG4dRICHSteppingAction() 
{
  // if the last hit was a zero energy deposit hit, it is
  // just reset and the memory is still allocated, so we
  // need to delete it here if the last hit was saved, hit
  // is a nullptr pointer which are legal to delete (it
  // results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
// This is the implementation of the G4 UserSteppingAction
bool EICG4dRICHSteppingAction::UserSteppingAction(const G4Step *aStep,
                                             bool was_used) 
{

  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
  {
    std::cout << "[>>>>>] call EICG4dRICHSteppingAction::UserSteppingAction" << std::endl;
  }

  // get touchables, points, volumes
  G4TouchableHandle preTouch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle postTouch = aStep->GetPostStepPoint()->GetTouchableHandle();
  G4VPhysicalVolume *preVol = preTouch->GetVolume();
  G4VPhysicalVolume *postVol = postTouch->GetVolume();
  G4StepPoint *prePoint = aStep->GetPreStepPoint();
  G4StepPoint *postPoint = aStep->GetPostStepPoint();

  // skip this step, if leaving the world (postVol will be nullptr)
  if (postPoint->GetStepStatus() == fWorldBoundary) 
  {
    if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "... skip this step (leaving world)" << std::endl;
    if (m_Hit) m_Hit->Reset();
    return false;
  }

  // get volume names
  G4String prePointVolName = prePoint->GetPhysicalVolume()->GetName();
  G4String postPointVolName = postPoint->GetPhysicalVolume()->GetName();
  G4String preTouchVolName = preVol->GetName();
  G4String postTouchVolName = postVol->GetName();

  // get track
  const G4Track *aTrack = aStep->GetTrack();
  G4String particleName = aTrack->GetParticleDefinition()->GetParticleName();
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
  {
    std::cout << "[-] track ID=" << aTrack->GetTrackID()
              << ", particle=" << particleName << std::endl;
  }

  // IsInDetector(preVol) returns
  //  == 0 outside of detector
  //  > 0 for hits in active volume
  //  < 0 for hits in passive material
  int whichactive = m_Detector->IsInDetector(preVol);
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
  {
    std::cout << "[_] step preVol=" << preTouchVolName
              << ", postVol=" << postTouchVolName << ", whichactive=" << whichactive
              << std::endl;
  }

  // reset hit classifiers
  hitType = -1;
  hitSubtype = -1;

  // classify hit type
  if (prePointVolName.contains("DRICHpetal") && postPointVolName.contains("DRICHpsst"))
  {
    hitType = hPSST;
  }
  else if (prePointVolName.contains("World") && postPointVolName.contains("DRICHvessel"))
  {
    hitType = hEntrance;
  }
  else if (prePointVolName.contains("DRICHvessel") && postPointVolName.contains("World"))
  {
    hitType = hExit;
  }
  else
  {
    hitType = hIgnore;
  }

  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE && hitType == hEntrance) 
  {  
    std::cout << "[__] step is ENTERING vessel" << std::endl;
  }

  // skip this step, if it's outside the detector, and not an entrance
  // or exit of the vessel
  if (!whichactive && hitType != hEntrance && hitType != hExit) 
  {
    if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "... skip this step" << std::endl;
    return false;
  }

  // get step energy // TODO: do we need `eion`?
  G4double edep = 0;
  G4double eion = 0;
  if (hitType != hEntrance) 
  {
    edep = aStep->GetTotalEnergyDeposit() / GeV;
    eion = (aStep->GetTotalEnergyDeposit() -
            aStep->GetNonIonizingEnergyDeposit()) / GeV;
  }
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[_] step edep=" << edep << ",   eion=" << eion << std::endl;

  // Here we have to decide if we need to create a new hit.  Normally this
  // should only be neccessary if a G4 Track enters a new volume or is freshly
  // created For this we look at the step status of the prePoint (beginning of
  // the G4 Step). This should be either fGeomBoundary (G4 Track crosses into
  // volume) or fUndefined (G4 Track newly created) Sadly over the years with
  // different G4 versions we have observed cases where G4 produces "impossible
  // hits" which we try to catch here These errors were always rare and it is
  // not clear if they still exist but we still check for them for safety. We
  // can reproduce G4 runs identically (if given the sequence of random number
  // seeds you find in the log), the printouts help us giving the G4 support
  // information about those failures
  switch (prePoint->GetStepStatus()) {

  // --- abnormal cases
  case fPostStepDoItProc: // step defined by PostStepDoItVector
    if (m_SavePostStepStatus != fGeomBoundary) 
    {
      // this is the okay case, fPostStepDoItProc called in a volume, not first
      // thing inside a new volume, just proceed here
      if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
        std::cout << "[__] first step in a new volume" << std::endl;
    } 
    else 
    {
      if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[ + ] step was defined by PostStepDoItVector" << std::endl;
      if (hitType != hEntrance) 
      {
        // this is an impossible G4 Step print out diagnostic to help debug, not
        // sure if this is still with us
        std::cerr << "ERROR: impossible G4 Step" << std::endl;
        std::cout << GetName() << ": New Hit for  " << std::endl;
        std::cout << "prestep status: "
                  << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
                  << ", poststep status: "
                  << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
                  << ", last pre step status: "
                  << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
                  << ", last post step status: "
                  << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus)
                  << std::endl;
        std::cout << "last track: " << m_SaveTrackId
                  << ", current trackid: " << aTrack->GetTrackID() << std::endl;
        std::cout << "phys pre vol: " << preTouchVolName
                  << " post vol : " << postTouch->GetVolume()->GetName() << std::endl;
        std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
                  << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
      } 
      else 
      {
        hitSubtype = entPostStep;
      }
    }

    // if this step is incident on the vessel, and we have not yet created a
    // hit, create one
    if (hitType == hEntrance) 
    {
      m_Hit = nullptr; // kill any leftover hit
      if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[++++] NEW hit (entrance)" << std::endl;
      m_Hit = new EICG4dRICHHit();
      this->InitHit(prePoint, aTrack, true);
    }
    break;

  // --- normal cases
  case fGeomBoundary:
  case fUndefined:
  default:

    // do nothing if not geometry boundary, not undefined, and not entrance
    if (prePoint->GetStepStatus() != fGeomBoundary &&
        prePoint->GetStepStatus() != fUndefined && hitType != hEntrance) 
    {
      if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[+] prepoint status ignored" << std::endl;
      break;
    }

    // create new hit
    if (!m_Hit) 
    {
      if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[++++] NEW hit" << std::endl;
      m_Hit = new EICG4dRICHHit();
      this->InitHit(prePoint, aTrack, true);
    } 
    else 
    {
      // if hit already exists, then Reset() has likely just been
      // called; in this case, initialize, but don't reset accumulators
      this->InitHit(prePoint, aTrack, false);
    }

    // print info
    if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) 
    {
      std::cout << "[+] prepoint status=";
      if (prePoint->GetStepStatus() == fGeomBoundary) std::cout << "fGeomBoundary";
      else if (prePoint->GetStepStatus() == fUndefined) std::cout << "fUndefined";
      else std::cout << "UNKNOWN!";
      std::cout << std::endl;
      if (aTrack->GetTrackID() > 1) 
      {
        std::cout << "[-] secondary track, creator process="
             << aTrack->GetCreatorProcess()->GetProcessName();
      } 
      else 
      {
        std::cout << "[-] primary track, particle=" << particleName;
      }
      std::cout << std::endl;
    }

    // tracking of the truth info // TODO: not used yet?
    if (G4VUserTrackInformation *p = aTrack->GetUserInformation()) 
    {
      if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p)) 
      {
        m_Hit->set_trkid(pp->GetUserTrackId());
        pp->GetShower()->add_g4hit_id(m_SaveHitContainer->GetID(), m_Hit->get_hit_id());
      }
    }
    break;
  }

  // This section is called for every step -------------------------------------
  // - some sanity checks for inconsistencies (aka bugs) we have seen over the years
  // - check if this hit was created, if not print out last post step status
  if (!m_Hit || !std::isfinite(m_Hit->get_x(0))) 
  {
    std::cout << "m_Hit @" << static_cast<void*>(m_Hit) << " isfinite = " << std::isfinite(m_Hit->get_x(0)) << std::endl;
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
    std::cout << "phys pre vol: " << preTouchVolName
              << " post vol : " << postTouch->GetVolume()->GetName() << std::endl;
    std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
              << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
    // This is fatal - a hit from nowhere. This needs to be looked at and fixed
    gSystem->Exit(1);
  }
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
  m_SaveVolPre = preVol;
  m_SaveVolPost = postVol;

  // update accumulators
  m_EdepSum += edep;
  if (whichactive > 0) 
  {
    m_EionSum += eion;
  }
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) 
  {
    std::cout << "[_] m_EdepSum=" << m_EdepSum << ",   m_EionSum=" << m_EionSum << std::endl;
  }

  // get petal number // TODO: does not work for entrance/exit hits
  // since they are identified by world<->vessel
  // crossings; the vessel has no petal number;
  // note that there is currently a gap between
  // the vessel and the petal volumes, which may
  // be unrealistic
  int petal = hitType == hEntrance ? m_Detector->GetPetal(postVol)
                                   : m_Detector->GetPetal(preVol);

  // save the hit ---------------------------------------------
  // if any of these conditions is true, this is the last step in
  // this volume and we consider saving the hit
  if (postPoint->GetStepStatus() == fGeomBoundary      /*left volume*/
      || postPoint->GetStepStatus() == fWorldBoundary  /*left world*/
      || postPoint->GetStepStatus() == fAtRestDoItProc /*track stops*/
      || aTrack->GetTrackStatus() == fStopAndKill)      /*track ends*/
  {

    if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
    {
      std::cout << "[---+] last step in the volume (pre=" << prePointVolName << ", post=" << postPointVolName << ")" << std::endl;
    }

    // hits to keep +++++++++++++++++++++++
    if (hitType != hIgnore) 
    {

      if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[-+] " << hitTypeStr[hitType] << " hit, KEEP!" << std::endl;

      // classify hit subtype
      switch (hitType) 
      {
        case hEntrance:
          if (hitSubtype != entPostStep) 
          {
            if (aTrack->GetTrackID() == 1) hitSubtype = entPrimary;
            else hitSubtype = entSecondary;
          }
          break;
        case hExit:
          if (aTrack->GetTrackID() == 1) hitSubtype = exPrimary;
          else hitSubtype = exSecondary;
          break;
        case hPSST:
          if (particleName == "opticalphoton") hitSubtype = psOptical;
          else if (particleName == "gamma") hitSubtype = psGamma;
          else hitSubtype = psOther;
          break;
        default:
          hitSubtype = subtypeUnknown;
      }
      if (hitSubtype == -1) hitSubtype = subtypeUnknown;

      // set hit vars
      m_Hit->set_hit_type_name(hitTypeStr[hitType]);
      m_Hit->set_hit_subtype_name(hitSubtypeStr[hitSubtype]);
      m_Hit->set_petal(petal);
      m_Hit->set_psst(m_Detector->GetPSST(postVol));
      m_Hit->set_pdg(aTrack->GetParticleDefinition()->GetPDGEncoding());
      m_Hit->set_particle_name(particleName);

      switch (hitType) 
      {
        case hEntrance:
          if (hitSubtype == entPostStep) m_Hit->set_process("postStep");
          else if (hitSubtype == entPrimary) m_Hit->set_process("primary");
          else  m_Hit->set_process(aTrack->GetCreatorProcess()->GetProcessName());
          break;
        case hExit:
          m_Hit->set_process("exitProcess");
          break;
        default:
          m_Hit->set_process(aTrack->GetCreatorProcess()->GetProcessName());
      }

      m_Hit->set_parent_id(aTrack->GetParentID());
      m_Hit->set_position(1, postPoint->GetPosition() / cm);
      m_Hit->set_momentum(aTrack->GetMomentum() / GeV);
      m_Hit->set_momentum_dir(aTrack->GetMomentumDirection() /*unitless*/);
      m_Hit->set_vertex_position(aTrack->GetVertexPosition() / cm);
      m_Hit->set_vertex_momentum_dir(aTrack->GetVertexMomentumDirection() /*unitless*/);
      m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

      // tracking of the truth info // TODO: not used yet?
      if (G4VUserTrackInformation *p = aTrack->GetUserInformation()) 
      {
        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p)) 
        {
          pp->SetKeep(1); // we want to keep the track
        }
      }

      // total accumulators
      m_Hit->set_edep(m_EdepSum);
      m_Hit->set_eion(m_EionSum);
      if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
      {
        std::cout << "[-+] m_EdepSum=" << m_EdepSum << ",    m_EionSum=" << m_EionSum << std::endl;
      }

      // transfer ownership to container
      m_SaveHitContainer->AddHit(petal, m_Hit);
      m_Hit = nullptr; // so that next track will create new hit
    }
    else 
    {
      // do not save this hit ++++++++++++++++++++++++++++
      // - reset hit object for reuse
      // - if last hit overall, memory still allocated
      // - does not reset local vars, such as m_EdepSum
      if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[-+] not keeping this hit" << std::endl;
      m_Hit->Reset();
    }
  }

  // return true to indicate the hit was used
  return true;
}

//____________________________________________________________________________..
void EICG4dRICHSteppingAction::InitHit(const G4StepPoint *prePoint_,
                                  const G4Track *aTrack_,
                                  bool resetAccumulators) 
{
  // set some entrance attributes, and track ID
  m_Hit->set_position(0, prePoint_->GetPosition() / cm);
  m_Hit->set_t(0, prePoint_->GetGlobalTime() / nanosecond);
  m_Hit->set_trkid(aTrack_->GetTrackID());
  m_SaveTrackId = aTrack_->GetTrackID();
  m_SaveHitContainer = m_HitContainer;
  // initializate accumulators (e.g., for total energy deposition)
  if (resetAccumulators) 
  {
    m_EdepSum = 0;
    m_EionSum = 0;
  }
}

//____________________________________________________________________________..
void EICG4dRICHSteppingAction::SetInterfacePointers(PHCompositeNode *topNode) 
{
  std::string hitnodename = "G4HIT_" + m_Detector->GetName();
  // now look for the map and grab a pointer to it.
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  // if we do not find the node we need to make it.
  if (!m_HitContainer) 
  {
    std::cout << "EICG4dRICHSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
}
