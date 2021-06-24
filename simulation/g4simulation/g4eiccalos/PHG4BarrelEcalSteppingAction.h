// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PHG4BARRELECALSTEPPINGACTION_H
#define G4DETECTORS_PHG4PHG4BARRELECALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4BarrelEcalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4BarrelEcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4BarrelEcalSteppingAction(PHG4BarrelEcalDetector*, const PHParameters* parameters);

  //! destroctor
  virtual ~PHG4BarrelEcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:

PHG4BarrelEcalDetector* m_Detector = nullptr;

  //! pointer to hit container
  PHG4HitContainer* m_HitContainer = nullptr;
  PHG4HitContainer* m_AbsorberHitContainer = nullptr;
  PHG4HitContainer* m_SupportHitContainer = nullptr;
  PHG4HitContainer* m_SaveHitContainer = nullptr;
  PHG4Hit* m_Hit = nullptr;
  PHG4Shower* m_SaveShower = nullptr;

  int m_ActiveFlag = 0;
  int m_AbsorberTruthFlag = 0;
  int m_SupportTruthFlag = 0;
  int m_BlackHoleFlag = 0;

  std::string m_HitNodeName;
  std::string m_AbsorberNodeName;
  std::string m_SupportNodeName;
  
};

#endif  // G4DETECTORS_PHG4BarrelEcalSTEPPINGACTION_H
