// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4RPSTEPPINGACTION_H
#define EICG4RPSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>
#include <string>

class EICG4RPDetector;
class EICG4RPSubsystem;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4Shower;
class PHG4HitContainer;
class PHParameters;

class EICG4RPSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  EICG4RPSteppingAction(EICG4RPSubsystem* subsys, EICG4RPDetector* detector, const PHParameters* parameters);

  //! destructor
  virtual ~EICG4RPSteppingAction() override;

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*) override;

  virtual void SaveLightYield(const int i = 1) { m_SaveLightYieldFlag = i; }

  virtual bool hasMotherSubsystem() const;

  virtual void SaveAllHits(bool i = true) { m_SaveAllHitsFlag = i; }

  virtual void HitNodeName(const std::string& name) { m_HitNodeName = name; }
  virtual void HitNodeNameVirt(const std::string& name) { m_HitNodeNameVirt = name; }

 private:
  //! pointer to the detector
  EICG4RPSubsystem* m_Subsystem;
  EICG4RPDetector* m_Detector;

  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
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
  std::string m_HitNodeName;
  std::string m_HitNodeNameVirt;
};

#endif  // EICG4RPSTEPPINGACTION_H
