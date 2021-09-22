// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BarrelEcalDETECTOR_H
#define G4DETECTORS_PHG4BarrelEcalDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Types.hh>  // for G4double
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Trap.hh>

#include <map>
#include <set>
#include <string>

class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4BarrelEcalDisplayAction;
class PHG4Subsystem;
class PHParameters;
class PHG4GDMLConfig;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4BarrelEcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4BarrelEcalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4BarrelEcalDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //!@name volume accessors
  int IsInBarrelEcal(G4VPhysicalVolume *) const;

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  int get_Layer() const { return m_Layer; }

 private:

  //! BECAL parameters

  G4double Radius = 85.*cm; //Inner radius of BECAL
  G4double tower_length = 45.5*cm; // Length of the Tower

  const int nTowers_layer = 128.; //Number of towers per phi tower
  const double becal_length = 415*cm;  //support width
  const double th = 1.;
  const double overlap = 0.1;

  G4double silicon_width_half;
  G4double kapton_width_half;
  G4double SIO2_width_half;
  G4double Carbon_width_half;
  G4double support_length;
  

  
  int PlaceTower(G4LogicalVolume *envelope);
  int ParseParametersFromTable();

  struct towerposition
  {

    G4double sizex1;
    G4double sizey1;
    G4double sizex2;
    G4double sizey2;
    G4double sizez;
    G4double pTheta;
    G4double centerx;
    G4double centery;
    G4double centerz;
    G4double rotx;
    G4double roty;
    G4double rotz;
    int idx_j;
    int idx_k;
  };

  G4Material *GetCarbonFiber();
  G4Material *GetSciGlass();

  G4Trap *GetTowerTrap(std::map<std::string, towerposition>::iterator iterator);
  G4Trap *GetSiTrap(std::map<std::string, towerposition>::iterator iterator);
  G4Trap* GetGlassTrap(std::map<std::string, towerposition>::iterator iterator);
  G4Trap* GetKaptonTrap(std::map<std::string, towerposition>::iterator iterator);
  G4Trap* GetSIO2Trap(std::map<std::string, towerposition>::iterator iterator);
  G4Trap* GetCarbonTrap(std::map<std::string, towerposition>::iterator iterator);

  G4LogicalVolume *GetTowerSci(std::map<std::string, towerposition>::iterator iterator);
  G4Trap *GetGlassTrapSubtract(std::map<std::string, towerposition>::iterator iterator);
  G4LogicalVolume *ConstructTower(std::map<std::string, towerposition>::iterator iterator);
  G4LogicalVolume *ConstructGlass(std::map<std::string, towerposition>::iterator iterator);
  G4LogicalVolume *ConstructSi(std::map<std::string, towerposition>::iterator iterator);
  G4LogicalVolume *ConstructKapton(std::map<std::string, towerposition>::iterator iterator);
  G4LogicalVolume *ConstructSIO2(std::map<std::string, towerposition>::iterator iterator);  
  G4LogicalVolume *ConstructCarbon(std::map<std::string, towerposition>::iterator iterator);  

  PHG4BarrelEcalDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;

  int m_ActiveFlag = 1;
  int m_AbsorberActiveFlag = 0;
  int m_SupportActiveFlag = 0;
  int m_Layer = 0;

 
  //G4LogicalVolume* singletower; 
  std::string m_TowerLogicNamePrefix;
  std::string m_SuperDetector;

  std::map<std::string, G4double> m_GlobalParameterMap;
  std::map<std::string, towerposition> m_TowerPostionMap;

  std::set<G4LogicalVolume *> m_AbsorberLogicalVolSet;
  std::set<G4LogicalVolume *> m_ScintiLogicalVolSet;
  std::set<G4LogicalVolume *> m_SupportLogicalVolSet;

  //! registry for volumes that should not be exported, i.e. fibers
  PHG4GDMLConfig* gdml_config = nullptr;
};

#endif
