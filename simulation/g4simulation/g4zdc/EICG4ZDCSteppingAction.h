// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4ZDCSTEPPINGACTION_H
#define EICG4ZDCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class EICG4ZDCDetector;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;
class PHG4Shower;

class EICG4ZDCSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  EICG4ZDCSteppingAction(EICG4ZDCDetector*, const PHParameters* parameters);

  //! destructor
  virtual ~EICG4ZDCSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  EICG4ZDCDetector* m_Detector;
  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
  PHG4Hit* m_Hit;
  PHG4HitContainer* m_SaveHitContainer;
  G4VPhysicalVolume* m_SaveVolPre;
  G4VPhysicalVolume* m_SaveVolPost;
  PHG4Shower* m_SaveShower;

  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_ActiveFlag;
  int m_BlackHoleFlag;
  double m_EdepSum;
  double m_EionSum;
  double m_LightYield;
 
};

#endif // EICG4ZDCSTEPPINGACTION_H
