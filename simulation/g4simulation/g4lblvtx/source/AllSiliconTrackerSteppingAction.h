// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ALLSILICONTRACKERSTEPPINGACTION_H
#define ALLSILICONTRACKERSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class AllSiliconTrackerDetector;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class AllSiliconTrackerSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  AllSiliconTrackerSteppingAction(AllSiliconTrackerDetector*, const PHParameters* parameters);

  //! destructor
  virtual ~AllSiliconTrackerSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

 private:
  //! pointer to the detector
  AllSiliconTrackerDetector* m_Detector;
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
  int m_BlackHoleFlag;
  double m_EdepSum;
};

#endif  // ALLSILICONTRACKERSTEPPINGACTION_H
