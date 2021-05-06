// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4EICDIRCSTEPPINGACTION_H
#define G4EICDIRCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class G4EicDircDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class G4EicDircSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  G4EicDircSteppingAction(G4EicDircDetector*, const PHParameters* parameters);

  //! destructor
  virtual ~G4EicDircSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  G4EicDircDetector* m_Detector;
  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
  PHG4Hit* m_Hit;
  PHG4HitContainer* m_SaveHitContainer;

  G4VPhysicalVolume* m_SaveVolPre;
  G4VPhysicalVolume* m_SaveVolPost;
  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_ActiveFlag;
  int m_BlackHoleFlag;
  double m_EdepSum;
  double m_EionSum;
};

#endif  // G4EICDIRCSTEPPINGACTION_H
