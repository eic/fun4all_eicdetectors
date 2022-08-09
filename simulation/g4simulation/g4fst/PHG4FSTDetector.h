// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FSTDETECTOR_H
#define G4DETECTORS_PHG4FSTDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4Types.hh>   // for G4double
#include <Geant4/G4Material.hh>

#include <map>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4FSTDisplayAction;
class PHG4FSTSteppingAction;
class PHG4Subsystem;
class PHParameters;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4FSTDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4FSTDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4FSTDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //!@name volume accessors
  int IsInActiveSensorFST(G4VPhysicalVolume *) const;

  //! Select mapping file for calorimeter tower
  void SetTowerMappingFile(const std::string &filename)
  {
    _mapping_tower_file = filename;
  }

  void SetTowerDimensions(G4double dx, G4double dy, G4double dz)
  {
    _tower_dx = dx;
    _tower_dy = dy;
    _tower_dz = dz;
  }

  void SetPlace(G4double place_in_x, G4double place_in_y, G4double place_in_z)
  {
    _place_in_x = place_in_x;
    _place_in_y = place_in_y;
    _place_in_z = place_in_z;
  }

  void SetXRot(G4double rot_in_x) { _rot_in_x = rot_in_x; }
  void SetYRot(G4double rot_in_y) { _rot_in_y = rot_in_y; }
  void SetZRot(G4double rot_in_z) { _rot_in_z = rot_in_z; }

  void SetActive(const int i = 1) { _active = i; }
  void SetAbsorberActive(const int i = 1) { _absorberactive = i; }

  int IsActive() const { return _active; }

  void SuperDetector(const std::string &name) { _superdetector = name; }
  void SetSteppingAction(PHG4FSTSteppingAction *stpact) { m_SteppingAction = stpact; }
  const std::string SuperDetector() const { return _superdetector; }

  int get_Layer() const { return _layer; }

  void BlackHole(const int i = 1) { _blackhole = i; }
  int IsBlackHole() const { return _blackhole; }

 private:
  void ConstructSTDisk(G4LogicalVolume *mother);
  void ConstructFST(G4LogicalVolume *mother);
  G4Material *MakeCarbonFoamMaterial_Longeron();
  G4Material *MakeCarbonFleece();
  G4Material *GetKapton();
  G4Material *MakeGlue();
  G4Material *MakeCarbonFoamMaterial_Wheel();
  G4Material *GetCarbonFiber();
  G4Material *MakeCarbonHoneyCombMaterial();
  G4Material *MakeCarbonFoamMaterial();

  PHG4FSTDisplayAction *m_DisplayAction;
  PHG4FSTSteppingAction *m_SteppingAction;

  /* Calorimeter envelope geometry */
  G4double _place_in_x;
  G4double _place_in_y;
  G4double _place_in_z;
  G4double _center_offset_x;
  G4double _center_offset_y;
  int _quadratic_detector;

  G4double _rot_in_x;
  G4double _rot_in_y;
  G4double _rot_in_z;

  G4double _rMin1;
  G4double _rMax1;
  G4double _rMin2;
  G4double _rMax2;

  G4double _dZ;
  G4double _sPhi;
  G4double _dPhi;

  /* DRCALO tower geometry */
  int _tower_type;
  G4double _tower_readout;
  G4double _tower_dx;
  G4double _tower_dy;
  G4double _tower_dz;

  G4double _scintFiber_diam;
  G4double _cerenkovFiber_diam;
  int _cerenkovFiber_material;
  int _tower_makeNotched;
  int _absorber_Material;

  G4double _wls_dw;
  G4double _support_dw;

  G4String _materialScintillator;
  G4String _materialAbsorber;

  int _active;
  int _absorberactive;
  int _layer;
  int _blackhole;

  std::string _towerlogicnameprefix;
  std::string _superdetector;
  std::string _mapping_tower_file;
  std::map<std::string, G4double> m_GlobalParameterMap;


  PHParameters *m_Params = nullptr;

};

#endif
