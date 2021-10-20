// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4B0STEPPINGACTION_H
#define EICG4B0STEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>
#include <string>

class EICG4B0Detector;
class EICG4B0Subsystem;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4Shower;
class PHG4HitContainer;
class PHParameters;

class EICG4B0SteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  EICG4B0SteppingAction(EICG4B0Subsystem* subsys, EICG4B0Detector* detector, const PHParameters* parameters);

  //! destructor
  virtual ~EICG4B0SteppingAction() override;

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*) override;

  virtual void SaveLightYield(const int i = 1) { m_SaveLightYieldFlag = i; }

  virtual bool hasMotherSubsystem() const;

  virtual void SaveAllHits(bool i = true) { m_SaveAllHitsFlag = i; }

  virtual void HitNodeName(const std::string& name) { m_HitNodeName = name; }

 private:
  //! pointer to the detector
  EICG4B0Subsystem* m_Subsystem;
  EICG4B0Detector* m_Detector;

  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
  PHG4Hit* m_Hit;
  //  PHG4HitContainer* m_SaveHitContainer;
  PHG4Shower* m_SaveShower;
  G4VPhysicalVolume* m_SaveVolPre;
  G4VPhysicalVolume* m_SaveVolPost;

  bool m_SaveAllHitsFlag = true;
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
};

#endif  // EICG4B0STEPPINGACTION_H
