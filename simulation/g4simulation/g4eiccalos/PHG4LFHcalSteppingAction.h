// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4LFHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4LFHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4TouchableHandle.hh>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4LFHcalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4LFHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4LFHcalSteppingAction(PHG4LFHcalDetector*, const PHParameters* parameters);

  //! destructor
  virtual ~PHG4LFHcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  PHG4LFHcalDetector* m_Detector = nullptr;

  //! pointer to hit container
  PHG4HitContainer* m_HitContainer = nullptr;
  PHG4HitContainer* m_AbsorberHitContainer = nullptr;
  PHG4HitContainer* m_SaveHitContainer = nullptr;
  PHG4Hit* m_Hit = nullptr;
  PHG4Shower* m_SaveShower = nullptr;

  int m_ActiveFlag = 0;
  int m_AbsorberTruthFlag = 0;
  int m_BlackHoleFlag = 0;
  int m_NlayersPerTowerSeg = 10;
  double m_zdepthtailcatcher = 100000;
  int m_NLayerTailCatcherOffset = 0;
};

#endif  // G4DETECTORS_PHG4LFHCALSTEPPINGACTION_H
