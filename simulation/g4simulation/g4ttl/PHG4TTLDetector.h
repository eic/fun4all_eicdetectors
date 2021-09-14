// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TTLDETECTOR_H
#define G4DETECTORS_PHG4TTLDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4DisplacedSolid.hh>     // for G4DisplacedSolid
#include <Geant4/G4ExceptionSeverity.hh>  // for FatalException, JustWarning
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MaterialTable.hh>  // for G4MaterialTable
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>  // for pi
#include <Geant4/G4Sphere.hh>
#include <Geant4/G4String.hh>         // for G4String
#include <Geant4/G4SystemOfUnits.hh>  // for cm, um, perCent
#include <Geant4/G4ThreeVector.hh>    // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>    // for G4Transform3D, G4RotateX3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4int
#include <Geant4/globals.hh>  // for G4Exception

#include "PHG4TTLSteppingAction.h"

#include <map>
#include <set>
#include <utility>

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4TTLDisplayAction;
// class PHG4TTLSteppingAction;
class PHG4Subsystem;
class PHParameters;

class PHG4TTLDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4TTLDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  ~PHG4TTLDetector(void) override
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  //!@name volume accessors
  //@{
  bool IsInSectorActive(G4VPhysicalVolume *physvol);
  //@}

  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  void SetSteppingAction(PHG4TTLSteppingAction *stpact) { m_SteppingAction = stpact; }

  // void OverlapCheck(const bool chk = true) override
  // {
  //   PHG4Detector::OverlapCheck(chk);
  //   // PHG4SectorConstructor::OverlapCheck(chk);
  // }
  void
  OverlapCheck(bool check = true) override
  {
    overlapcheck_sector = check;
  }
  // void Verbosity(int v) override {m_Verbosity = v;}
  // int Verbosity() const {return m_Verbosity;}

 public:
  // properties

  std::string name_base;

 private:
  PHG4TTLDisplayAction *m_DisplayAction;
  PHG4TTLSteppingAction *m_SteppingAction;
  int m_Verbosity;
  std::string superdetector;
  PHParameters *m_Params = nullptr;

 protected:
  G4LogicalVolume *
  RegisterLogicalVolume(G4LogicalVolume *);
  bool overlapcheck_sector;

  G4VSolid *
  Construct_Sectors_Plane(           //
      const std::string &name,       //
      const double start_z,          //
      const double thickness,        //
      G4VSolid *SecConeBoundary_Det  //
  );
  typedef std::map<G4String, G4LogicalVolume *> map_log_vol_t;
  map_log_vol_t map_log_vol;

  G4PVPlacement *
  RegisterPhysicalVolume(G4PVPlacement *v, const bool active = false);

  void BuildForwardTTL(G4LogicalVolume *world);
  void BuildBarrelTTL(G4LogicalVolume *world);

  typedef std::pair<G4String, G4int> phy_vol_idx_t;
  typedef std::map<phy_vol_idx_t, G4PVPlacement *> map_phy_vol_t;
  map_phy_vol_t map_phy_vol;         //! all physics volume
  map_phy_vol_t map_active_phy_vol;  //! active physics volume
};

#endif
