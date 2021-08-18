#include "PrtOpBoundaryProcess.h"

#include "G4EicDircDetector.h"
#include <Geant4/G4ios.hh>
//#include <Geant4/G4TouchableHandle.hh>

PrtOpBoundaryProcess::PrtOpBoundaryProcess()
  : G4OpBoundaryProcess()
{}

G4VParticleChange* PrtOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{  
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  //G4TouchableHandle touch = pPreStepPoint->GetTouchableHandle();
  //G4TouchableHandle touchpost = pPostStepPoint->GetTouchableHandle();

  // get volume of the current step
  //G4VPhysicalVolume *volume = pPreStepPoint->GetPhysicalVolume();
  //G4VPhysicalVolume *volume_post = pPostStepPoint->GetPhysicalVolume();

  /*G4String vol_name = volume->GetName();

  G4cout << "vol name = " << vol_name << G4endl;

  int whichactive_int = m_Detector->IsInDetector(volume);
  //int whichactive_int_post = m_Detector->IsInDetector(volume_post);
  //int whichactive_int = G4EicDircDetector::IsInDetector(volume);

  //bool whichactive = (whichactive_int == 3);
  */
  G4VParticleChange* particleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep); 

  //if(whichactive_int == 1)  
    
  // int parentId = aTrack.GetParentID();
  // std::cout<<"parentId   "<<parentId <<std::endl;
  // if(parentId==1) particleChange->ProposeTrackStatus(fStopAndKill);

  //double endofbar = 0.5*(4200+4*0.05); //1250/2.;
  
  // LUT
  /*if(PrtManager::Instance()->GetRunType() == 1 && pPostStepPoint->GetPosition().z() > pPreStepPoint->GetPosition().z()){
    if(PrtManager::Instance()->GetEvType() != 1 ) particleChange->ProposeTrackStatus(fStopAndKill);
    if(pPreStepPoint->GetPosition().z() > endofbar) particleChange->ProposeTrackStatus(fStopAndKill);
    }*/

  if(aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wExpVol" && pPostStepPoint->GetPosition().z() > pPreStepPoint->GetPosition().z()){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3" && pPostStepPoint->GetPosition().z() > pPreStepPoint->GetPosition().z()){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
  
  
  // kill photons outside bar and prizm
  if(GetStatus() == FresnelRefraction 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens1" 
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens2")
     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // // black edge of the lens3
  // if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
  //     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc")
  //    || (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
  // 	 &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wLens3")){
  //   particleChange->ProposeTrackStatus(fStopAndKill);
  // }
  
  
  if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens1" 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wLens1"){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
  if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens2" 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wLens2"){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
    

  return particleChange;
    
}
