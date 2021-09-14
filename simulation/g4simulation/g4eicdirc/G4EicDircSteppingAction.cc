#include "G4EicDircSteppingAction.h"
#include "G4EicDircDetector.h"
#include "PrtHit.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4StepStatusDecode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <TSystem.h>
#include <TVector3.h>

#include <Geant4/G4NavigationHistory.hh>
#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>            // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>        // for G4TouchableHandle
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4TrackStatus.hh>            // for fStopAndKill
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/Randomize.hh>

#include <cmath>  // for isfinite
#include <iostream>
#include <string>  // for operator<<, string

class PHCompositeNode;
//____________________________________________________________________________..
G4EicDircSteppingAction::G4EicDircSteppingAction(
    G4EicDircDetector *detector, const PHParameters *parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_BlackHoleFlag(m_Params->get_int_param("blackhole"))
{
}

G4EicDircSteppingAction::~G4EicDircSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}


//____________________________________________________________________________..
bool G4EicDircSteppingAction::UserSteppingAction(const G4Step *aStep,
                                                 bool was_used)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume *volume = touch->GetVolume();
  G4VPhysicalVolume *volume_post = touchpost->GetVolume();
  // IsInDetector(volume) returns
  //  == 0 outside of detector
  //   > 0 for hits in active volume
  //  < 0 for hits in passive material
  int whichactive_int = m_Detector->IsInDetector(volume);
  int whichactive_int_post = m_Detector->IsInDetector(volume_post);
  bool whichactive = (whichactive_int > 0 && whichactive_int < 12);
  //int whichactive = m_Detector->IsInDetector(volume);
  if (!whichactive)
  {
    return false;
  }


  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion =
      (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) /
      GeV;
  const G4Track *aTrack = aStep->GetTrack();


  /*if(aTrack->GetCurrentStepNumber()>50000 || aTrack->GetTrackLength() > 30000) 
    {
      G4Track *killtrack0 = const_cast<G4Track *>(aTrack);
      killtrack0->SetTrackStatus(fStopAndKill);
      //return false;
      }*/


  /*if((whichactive_int == 10 && whichactive_int_post == 1) || (whichactive_int == 2 && whichactive_int_post == 1))
    {
      G4Track *killtrack = const_cast<G4Track *>(aTrack);
      killtrack->SetTrackStatus(fStopAndKill);
    }
  */
 
  // if this block stops everything, just put all kinetic energy into edep
  if (m_BlackHoleFlag)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track *killtrack = const_cast<G4Track *>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }


  int detector_id = 0;  // we use here only one detector in this simple example
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
  // cout << "track id " << aTrack->GetTrackID() << endl;
  //       cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
  //       cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;

  G4String prePointVolName = prePoint->GetPhysicalVolume()->GetName();
  G4String postPointVolName = postPoint->GetPhysicalVolume()->GetName();    

  switch (prePoint->GetStepStatus())
  {
  case fPostStepDoItProc:
    if (m_SavePostStepStatus != fGeomBoundary)
    {
      break;
    }
    else
    {
      // this is bad from G4 print out diagnostic to help debug, not sure if
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
 
  case fGeomBoundary:
  case fUndefined:
    if (!m_Hit)
      {
      //m_Hit = new PHG4Hitv1();
      m_Hit = new PrtHit();
      }
  
    // for momentum direction at bar
    //if((prePointVolName.contains("wBar")) && (aStep->IsFirstStepInVolume()) && (aTrack->GetParentID() == 0))
    
	//m_SaveHitContainer->AddHit(detector_id, m_Hit);
	// ownership has been transferred to container, set to null                                                                                  
    	// so we will create a new hit for the next track                                                                                            
     
    m_Hit->set_layer(detector_id);
    // here we set the entrance values in cm
    m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
    m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
    m_Hit->set_z(0, prePoint->GetPosition().z() / cm);
    // time in ns
    m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
    // set the track ID
    m_Hit->set_trkid(aTrack->GetTrackID());
    m_SaveTrackId = aTrack->GetTrackID();

    // set the initial energy deposit
    m_EdepSum = 0;
    if (whichactive > 0)
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
    // this is for the tracking of the truth info
    if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
      {
        m_Hit->set_trkid(pp->GetUserTrackId());
        pp->GetShower()->add_g4hit_id(m_SaveHitContainer->GetID(),
                                      m_Hit->get_hit_id());
      }
    }

    break;
  default:
    break;
  }
  //cout << "detector id = " << detector_id << endl;

  // some sanity checks for inconsistencies (aka bugs)
  // check if this hit was created, if not print out last post step status
  if (!m_Hit || !isfinite(m_Hit->get_x(0)))
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

    gSystem->Exit(1);
  }
  m_SavePreStepStatus = prePoint->GetStepStatus();
  m_SavePostStepStatus = postPoint->GetStepStatus();
  m_SaveVolPre = volume;
  m_SaveVolPost = touchpost->GetVolume();

  // here we just update the exit values, it will be overwritten
  // for every step until we leave the volume or the particle
  // ceases to exist
  // sum up the energy to get total deposited
  m_EdepSum += edep;
  if (whichactive > 0)
  {
    m_EionSum += eion;
  }

	    	        
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
      aTrack->GetTrackStatus() == fStopAndKill)
  {       
    //if((prePoint->GetStepStatus() == fGeomBoundary) && 
    if(whichactive_int == 9 || whichactive_int == 7 || whichactive_int == 8) // for relection information (7-wLens2, 8-wLens3, 9-wPrizm) 
      	{	 
	  G4String vname = touch->GetVolume()->GetName();
	     
	  // normal to the closest boundary
	  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

	  Int_t nid = 0;	
	  G4bool valid;
	  G4ThreeVector normal = theNavigator->GetLocalExitNormal(&valid);
	  G4ThreeVector gnormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(-normal);
	  normal = touch->GetHistory()->GetTransform(1).TransformAxis(gnormal); // in lDirc

	  //Int_t prizm_hit_trackid = aTrack->GetTrackID();
	  Int_t prizm_hit_trackid = m_SaveTrackId;
  
	  if (valid)
	    {
	      if(vname=="wPrizm")
		{
		  if(normal.y()> 0.99) nid = 1; // right
		  if(normal.y()<-0.99) nid = 2; // left
		  if(normal.x()>0.99) nid = 3; // bottom
		  if(fabs(normal.x()+0.866025)<0.1 ) nid = 4;
		}
	      else if(vname=="wLens3")
		{
		  if(normal.y()> 0.99) nid = 5; // right
		  if(normal.y()<-0.99) nid = 6; // left
		  if(normal.x()>0.99) nid = 7; // bottom
		  if(normal.x()<-0.99) nid = 8; // top
		}
	      else if(vname=="wLens2")
		{
		  nid = 9;
		}
	     
	      if(nid > 0)
		{		 
		  vector_trackid.push_back(prizm_hit_trackid);
		  vector_nid.push_back(nid);
		}
	     
	    }
	}   
  

   // save only hits with energy deposit (or geantino)
    
    if (m_EdepSum > 0 || geantino)
    {
      // update values at exit coordinates and set keep flag
      // of track to keep
      m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
      m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
      m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

      m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
            	      
      if(whichactive_int_post == 11) // post step in Pixel ---------------
	{
      // Get cell id 
      //G4int layerNumber = touchpost->GetReplicaNumber(0);
      //const G4DynamicParticle* dynParticle = aTrack->GetDynamicParticle();
      //G4ParticleDefinition* particle = dynParticle->GetDefinition();  
      //G4String ParticleName = particle->GetParticleName();
  
      G4ThreeVector globalpos = aStep->GetPostStepPoint()->GetPosition();
      G4ThreeVector localpos = touchpost->GetHistory()->GetTopTransform().TransformPoint(globalpos);
      G4ThreeVector translation = touchpost->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,0));
      G4ThreeVector inPrismpos = touchpost->GetHistory()->GetTransform( 1 ).TransformPoint(globalpos);
      G4ThreeVector g4mom = aTrack->GetVertexMomentumDirection();//GetMomentum();
      G4ThreeVector g4pos = aTrack->GetVertexPosition();
 
      //G4ThreeVector localvec = touchpost->GetHistory()->GetTopTransform().TransformAxis(g4mom);

      TVector3 globalPos(inPrismpos.x(),inPrismpos.y(),inPrismpos.z());
      TVector3 localPos(localpos.x(),localpos.y(),localpos.z());
      TVector3 digiPos(translation.x(),translation.y(),translation.z());
      TVector3 momentum(g4mom.x(),g4mom.y(),g4mom.z());
      TVector3 position(g4pos.x(),g4pos.y(),g4pos.z());

      // information from prizm
      double time = aStep->GetPreStepPoint()->GetLocalTime();

      int mcp = touchpost->GetReplicaNumber(1);
      int pix = touchpost->GetReplicaNumber(0);     
     
      Double_t wavelength = 1.2398/(aTrack->GetMomentum().mag()*1E6)*1000;
      
      // transport efficiency ----

      /*double pi(4 * atan(1));
      double roughness(0.5); // nm
      double angleX = localvec.angle(G4ThreeVector(1, 0, 0));
      double angleY = localvec.angle(G4ThreeVector(0, 1, 0));
      if (angleX > 0.5 * pi) angleX = pi - angleX;
      if (angleY > 0.5 * pi) angleY = pi - angleY;
      double length = aTrack->GetTrackLength() - 400; // 400 - average path in EV
      double lengthx = fabs(length * localvec.x());  // along the bar
      double lengthy = fabs(length * localvec.y());

      //cout << "track length = " << aTrack->GetTrackLength() << endl;

      int nBouncesX = (int)(lengthx) / 17;
      int nBouncesY = (int)(lengthy) / (358.5/11);

      double ll = wavelength * wavelength;
      double n_quartz = sqrt(1. + (0.696 * ll / (ll - pow(0.068, 2))) +
			     (0.407 * ll / (ll - pow(0.116, 2))) + 0.897 * ll / (ll - pow(9.896, 2)));
      double bounce_probX = 1 - pow(4 * pi * cos(angleX) * roughness * n_quartz / wavelength, 2);
      double bounce_probY = 1 - pow(4 * pi * cos(angleY) * roughness * n_quartz / wavelength, 2);

      double totalProb = pow(bounce_probX, nBouncesX) * pow(bounce_probY, nBouncesY);

      if (G4UniformRand() < totalProb) //{
      */     

      // time since track created
      m_Hit->SetLeadTime(time);
      m_Hit->SetTotTime(wavelength); //set photon wavelength
      m_Hit->SetMcpId(mcp);
      m_Hit->SetPixelId(pix);
      m_Hit->SetGlobalPos(globalPos);
      m_Hit->SetLocalPos(localPos);
      m_Hit->SetDigiPos(digiPos);
      m_Hit->SetPosition(position);
      m_Hit->SetMomentum(momentum);
                  
      int refl = 0;
      //Int_t normal_id = 0;
      Long64_t pathId = 0;
      //TVector3 mom_bar;
      //TVector3 pos_bar;

      for(std::vector<Int_t>::size_type i = 0; i < vector_trackid.size(); i++)
	{
	  //if(aTrack->GetTrackID() == vector_trackid[i]) 
	  if(m_SaveTrackId == vector_trackid[i])
	    {
	      ++refl;
	      Int_t normal_id = vector_nid[i];
		    //std::cout << "nid = " << normal_id << std::endl;
	      pathId = (pathId * 10L) + normal_id;
	    }
	}
	    
      //std::cout << "nrefl = " << refl << std::endl;
      //std::cout << "path id = " << pathId << std::endl;		  

      m_Hit->SetNreflectionsInPrizm(refl);
      m_Hit->SetPathInPrizm(pathId);
      
      /*for(std::vector<Int_t>::size_type i = 0; i < vector_bar_hit_trackid.size(); i++)
	{
	  if(aTrack->GetParentID() == vector_bar_hit_trackid[i])
	    {
	      mom_bar = vector_p_bar[i];
	      pos_bar = vector_hit_pos_bar[i];
	    }
	}
      
      m_Hit->SetMomentumAtBar(mom_bar);
      m_Hit->SetPositionAtBar(pos_bar);
      */

      //m_Hit->SetParticleId(aTrack->GetTrackID());
      //hit.SetParentParticleId(aTrack->GetParentID());
      

      
      //if((prePointVolName.contains("World")) && (postPointVolName.contains("wBar")) && (aTrack->GetParentID() == 0))
      //{
	/*if (!m_Hit)
	  {
	    m_Hit = new PrtHit();
	    }*/
	  
	/*G4ThreeVector momentum_at_bar = aTrack->GetMomentum();
	G4ThreeVector position_at_bar = prePoint->GetPosition();

	TVector3 p_bar(momentum_at_bar.x(), momentum_at_bar.y(), momentum_at_bar.z());
	TVector3 hit_pos_bar(position_at_bar.x(), position_at_bar.y(), position_at_bar.z());
	
	Int_t bar_hit_trackid = aTrack->GetTrackID();
	//detector_id = 1;
	
	//vector_bar_hit_trackid.push_back(bar_hit_trackid);
	//vector_p_bar.push_back(p_bar);
	//vector_hit_pos_bar.push_back(hit_pos_bar);

	//bar_vectors::vector_p_bar.push_back(p_bar); 
	//bar_vectors::vector_hit_pos_bar.push_back(hit_pos_bar);

	TVector3 mom_bar = p_bar;                                                                                                                  
	TVector3 pos_bar = hit_pos_bar;
	m_Hit->SetMomentumAtBar(mom_bar);                                                                                                          
      	m_Hit->SetPositionAtBar(pos_bar);
	*/
      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
      if (geantino)
      {
        m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way
                              // geantinos survive the g4hit compression
        if (whichactive > 0)
        {
          m_Hit->set_eion(-1);
        }
      }
      else
      {
        m_Hit->set_edep(m_EdepSum);
      }
      if (whichactive > 0)
      {
        m_Hit->set_eion(m_EionSum);
      }
      //} // pixel volume ends

      m_SaveHitContainer->AddHit(detector_id, m_Hit);
	    	
      // ownership has been transferred to container, set to null
      // so we will create a new hit for the next track
      //m_Hit = nullptr;
      m_Hit = nullptr;
	}
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
void G4EicDircSteppingAction::SetInterfacePointers(PHCompositeNode *topNode)
{
  if (!m_HitNodeName.empty())
  {
    m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  }
  if (!m_AbsorberNodeName.empty())
  {
    m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_AbsorberNodeName);
    if (!m_AbsorberHitContainer)
    {
      if (Verbosity() > 0)
      {
	std::cout << "G4EicDircSteppingAction::SetTopNode - unable to find " << m_AbsorberNodeName << std::endl;
      }
    }
  }
  if (! m_SupportNodeName.empty())
  {
    m_SupportHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_SupportNodeName);
    if (!m_SupportHitContainer)
    {
      if (Verbosity() > 0)
      {
	std::cout << "G4EicDircSteppingAction::SetTopNode - unable to find " << m_SupportNodeName << std::endl;
      }
    }

  }
  if (!m_HitContainer)
  {
    std::cout << "G4EicDircSteppingAction::SetTopNode - unable to find " << m_HitNodeName << std::endl;
  }

}




