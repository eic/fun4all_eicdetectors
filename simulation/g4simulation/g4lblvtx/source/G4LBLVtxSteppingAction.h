// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4LBLVTXSTEPPINGACTION_H
#define G4LBLVTXSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class G4LBLVtxDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class G4LBLVtxSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  G4LBLVtxSteppingAction(G4LBLVtxDetector*, const PHParameters* parameters);

  //! destructor
  virtual ~G4LBLVtxSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  G4LBLVtxDetector* m_Detector;

  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
  PHG4HitContainer* m_AbsorberHitContainer;
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

#endif  //  G4LBLVTXSTEPPINGACTION_H
