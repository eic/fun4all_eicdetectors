// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TRDSTEPPINGACTION_H
#define G4DETECTORS_PHG4TRDSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4TRDDetector;
class PHG4TRDSubsystem;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4TRDSteppingAction : public PHG4SteppingAction 
{

public:
  //! constructor
  // PHG4TRDSteppingAction(PHG4TRDSubsystem *subsys, PHG4TRDDetector *detector, const PHParameters *parameters);
   PHG4TRDSteppingAction(PHG4TRDDetector *detector, const PHParameters *parameters);

  //! destructor
  ~PHG4TRDSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  // needed for hit position crosschecks, if this volume is inside
  // another volume the absolut hit coordinates in our G4Hits and
  // the local coordinates differ, so checking against our place in z
  // goes wrong
  //bool hasMotherSubsystem() const;

  //void HitNodeName(const std::string &name) {m_HitNodeName = name;}
  
private:
  //! Pointer to subsystem
  PHG4TRDSubsystem *m_Subsystem;
  //! pointer to the detector
  PHG4TRDDetector *m_Detector;
  const PHParameters *m_Params;
  
  //! pointer to hit container
  PHG4HitContainer *m_HitContainer;
  PHG4HitContainer *m_ActiveGasHits;
  PHG4Hit *m_Hit;
  PHG4HitContainer *m_SaveHitContainer;
  PHG4Shower *m_SaveShower;
  
  G4VPhysicalVolume *m_SaveVolPre;
  G4VPhysicalVolume *m_SaveVolPost;
  
  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_BlackHoleFlag;
  int m_ActiveFlag;
  int m_UseG4StepsFlag;
  double m_Zmin;
  double m_Zmax;
  double m_EdepSum;
  //std::string m_HitNodeName;
  
  


};

#endif 
