// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4LUMISTEPPINGACTION_H
#define EICG4LUMISTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>
#include <string>

#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4StepPoint.hh> 


class EICG4LumiDetector;
class EICG4LumiSubsystem;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4Shower;
class PHG4HitContainer;
class PHParameters;

class EICG4LumiSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  EICG4LumiSteppingAction(EICG4LumiSubsystem *subsys, EICG4LumiDetector*, const PHParameters* parameters);

  //! destructor
  virtual ~EICG4LumiSteppingAction() override; 

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*) override;

  virtual void SaveLightYield(const int i = 1) { m_SaveLightYieldFlag = i;}
  virtual bool hasMotherSubsystem() const;
  virtual void SaveAllHits(bool i = true){ m_SaveAllHitsFlag = i;}
  virtual void HitNodeNameCAL(const std::string &name) {m_HitNodeNameCAL=name;}
  virtual void HitNodeNameTracking(const std::string &name) {m_HitNodeNameTracking=name;}
  virtual void HitNodeNameVirt(const std::string &name) {m_HitNodeNameVirt = name;}

 private:

  int FindTowerIndexFromPosition(G4StepPoint* prePoint, int& j, int& k);

  //! pointer to the detector
  EICG4LumiSubsystem* m_Subsystem;

  //! pointer to the detector
  EICG4LumiDetector* m_Detector;
  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainerCAL;
  PHG4HitContainer* m_HitContainerTracking;
  PHG4HitContainer* m_HitContainerVirt;
  PHG4Hit* m_Hit;
  PHG4Shower* m_SaveShower;
  G4VPhysicalVolume* m_SaveVolPre;
  G4VPhysicalVolume* m_SaveVolPost;

  bool m_SaveAllHitsFlag = false;
  int m_SaveLightYieldFlag;
  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_ActiveFlag; 
  int m_BlackHoleFlag; 
  int m_UseG4StepsFlag;
  double m_Zmin;
  double m_Zmax;
  double m_Tmin;
  double m_Tmax;
  double m_EdepSum;
  double m_EabsSum;
  double m_EionSum;
  std::string m_HitNodeNameCAL;
  std::string m_HitNodeNameTracking;
  std::string m_HitNodeNameVirt;

};

#endif // EICG4LUMISTEPPINGACTION_H
