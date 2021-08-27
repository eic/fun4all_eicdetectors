#include "G4EicDircOpBoundaryProcess.h"

#include <Geant4/G4Step.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4TouchableHistory.hh>

G4VParticleChange* G4EicDircOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{  
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  // int parentId = aTrack.GetParentID();
  // std::cout<<"parentId   "<<parentId <<std::endl;
  // if(parentId==1) pParticleChange->ProposeTrackStatus(fStopAndKill);

  double endofbar = 0.5*(4200+4*0.05); //1250/2.;
  
  // ideal focusing
//  if(PrtManager::Instance()->GetLens() == 10)
{
    G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
    G4TouchableHistory* touchable = (G4TouchableHistory*)(pPostStepPoint->GetTouchable());
    G4ThreeVector lpoint =  touchable->GetHistory()->GetTransform( 1 ).TransformPoint(theGlobalPoint1);
    if(lpoint.getZ() < endofbar+0.0001 && lpoint.getZ() > endofbar-0.0001){
      G4ThreeVector ww  = pPreStepPoint->GetTouchableHandle()->GetHistory()->
	GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,endofbar));

      if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()!="wGlue") 
	pParticleChange->ProposeTrackStatus(fStopAndKill);
      else
	aParticleChange.ProposePosition(ww.getX(), ww.getY(),lpoint.getZ()-0.0005);
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }

{  
 pParticleChange->ProposeTrackStatus(fStopAndKill);
    if(pPreStepPoint->GetPosition().z() < endofbar) pParticleChange->ProposeTrackStatus(fStopAndKill);
}

  if(aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wExpVol" && pPostStepPoint->GetPosition().z()<pPreStepPoint->GetPosition().z()){
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }

  if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3" && pPostStepPoint->GetPosition().z()<pPreStepPoint->GetPosition().z()){
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }
  
  
  // kill photons outside bar and prizm
/*
  if(GetStatus() == FresnelRefraction 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }
*/
  if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens1" 
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens2")
     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }

  // // black edge of the lens3
  // if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
  //     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc")
  //    || (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
  // 	 &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wLens3")){
  //   pParticleChange->ProposeTrackStatus(fStopAndKill);
  // }
  
  
  if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens1" 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wLens1"){
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }
  if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens2" 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wLens2"){
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }


  return pParticleChange;

}
