// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TTLSTEPPINGACTION_H
#define G4DETECTORS_PHG4TTLSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4ParticleDefinition.hh>  // for G4ParticleDefinition
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>              // for G4StepPoint
#include <Geant4/G4StepStatus.hh>             // for fGeomBoundary, fAtRestD...
#include <Geant4/G4String.hh>                 // for G4String
#include <Geant4/G4SystemOfUnits.hh>          // for cm, GeV, nanosecond
#include <Geant4/G4ThreeVector.hh>            // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>        // for G4TouchableHandle
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4TrackStatus.hh>            // for fStopAndKill
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation


#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>
#include <TVector3.h>

#include <iostream>
#include <string>  // for string, operator+, oper...

class G4Step;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHG4TTLDetector;
class PHG4Shower;

class PHG4TTLSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4TTLSteppingAction(PHG4TTLDetector*);

  //! destructor
  ~PHG4TTLSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;
  void SetLGADResolution(double lgadReso){
    _sensor_resolution_x = lgadReso;
    _sensor_resolution_y = lgadReso;
    };
  void SetNPhiModules(double nPhimod){
    _N_phi_modules = nPhimod;
    };
  void SetZPositionFwd(double zpos){
    _z_pos_TTL = zpos;
    };
  void SetIsForwardTTL(bool isfwd){
    _isFwd_TTL = isfwd;
    };
  void CalculateSensorHitIndices(G4StepPoint* prePoint, int& module_ret, int& layer, int& sensor_0, int& sensor_1, int& j, int& k, TVector3& sensorposition);

 private:
  //! pointer to the detector
  PHG4TTLDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4Hit* hit;
  PHG4Shower* saveshower;

  int layer_id;
  bool _isFwd_TTL;
  double _N_phi_modules;
  double _z_pos_TTL;
  double _sensor_resolution_x;
  double _sensor_resolution_y;
};

#endif  //__G4PHPHYTHIAREADER_H__
