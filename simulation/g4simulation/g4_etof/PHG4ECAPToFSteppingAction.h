#ifndef G4DETECTORS_PHG4ECAPToFSTEPPINGACTION_H
#define G4DETECTORS_PHG4ECAPToFSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ECAPToFDetector;
class PHG4ECAPToFSubsystem;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4ECAPToFSteppingAction : public PHG4SteppingAction 
{

public:
  //! constructor
  // PHG4ECAPToFSteppingAction(PHG4ECAPToFSubsystem *subsys, PHG4ECAPToFDetector *detector, const PHParameters *parameters);
   PHG4ECAPToFSteppingAction(PHG4ECAPToFDetector *detector, const PHParameters *parameters);

  //! destructor
  ~PHG4ECAPToFSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  // needed for hit position crosschecks, if this volume is inside
  // another volume the absolut hit coordinates in our G4Hits and
  // the local coordinates differ, so checking against our place in z
  // goes wrong
  //bool hasMotherSubsystem() const;


  /// this is just needed for use as reference plane for projections
  // this is the only detector using this - there is no need to add
  // this to our parameters
  void SaveAllHits(bool i = true) { m_SaveAllHitsFlag = i; }
  //void HitNodeName(const std::string &name) {m_HitNodeName = name;}
  
private:
  //! Pointer to subsystem
  PHG4ECAPToFSubsystem *m_Subsystem;
  //! pointer to the detector
  PHG4ECAPToFDetector *m_Detector;
  const PHParameters *m_Params;
  
  //! pointer to hit container
  PHG4HitContainer *m_HitContainer;
  //PHG4HitContainer *m_AbsorberHits;
  PHG4HitContainer *m_ActiveGasHits;
  PHG4Hit *m_Hit;
  PHG4HitContainer *m_SaveHitContainer;
  PHG4Shower *m_SaveShower;
  
  G4VPhysicalVolume *m_SaveVolPre;
  G4VPhysicalVolume *m_SaveVolPost;
  bool m_SaveAllHitsFlag = false;
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
