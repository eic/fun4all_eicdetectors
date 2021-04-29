// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDDUALREADOUTSTEPPINGACTION_H
#define G4DETECTORS_PHG4FORWARDDUALREADOUTSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4StepPoint.hh>               // for G4StepPoint

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ForwardDualReadoutDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4ForwardDualReadoutSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4ForwardDualReadoutSteppingAction(PHG4ForwardDualReadoutDetector*, const int absorberactive);

  //! destructor
  virtual ~PHG4ForwardDualReadoutSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);
  void SetTowerSize(G4double twrsze)
    {
      _tower_size = twrsze;
    }
  void SetDetectorSize(G4double detsze)
    {
      _detector_size = detsze;
    }
 private:
  int FindTowerIndex(G4TouchableHandle touch, int& j, int& k);
  int FindTowerIndexFromPosition(G4StepPoint* prePoint, int& j, int& k);

  int ParseG4VolumeName(G4VPhysicalVolume* volume, int& j, int& k);

  //! pointer to the detector
  PHG4ForwardDualReadoutDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4HitContainer* absorberhits_;
  PHG4HitContainer* hitcontainer;
  PHG4Hit* hit;
  PHG4Shower* saveshower;

  G4double _tower_size;
  G4double _detector_size;
  int absorbertruth;
  int light_scint_model;
};

#endif  // G4DETECTORS_PHG4FORWARDDUALREADOUTSTEPPINGACTION_H
