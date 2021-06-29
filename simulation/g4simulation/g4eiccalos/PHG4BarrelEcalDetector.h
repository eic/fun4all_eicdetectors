// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BarrelEcalDETECTOR_H
#define G4DETECTORS_PHG4BarrelEcalDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Types.hh>  // for G4double
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>

#include <map>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4BarrelEcalDisplayAction;
class PHG4Subsystem;
class PHParameters;

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

  G4LogicalVolume *ConstructTower();
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower);
  int ParseParametersFromTable();

  virtual std::pair<G4LogicalVolume*, G4Transform3D>
  Construct_AzimuthalSeg();

  struct towerposition
  {

    G4double centerx;
    G4double centery;
    G4double centerz;
    G4double roty;
    G4double rotz;
    int idx_j;
    int idx_k;
  };

  PHG4BarrelEcalDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;

  int m_ActiveFlag = 1;
  int m_AbsorberActiveFlag = 0;
  int m_SupportActiveFlag = 0;
  int m_Layer = 0;

 
   G4LogicalVolume* singletower; 
  std::string m_TowerLogicNamePrefix;
  std::string m_SuperDetector;

  std::map<std::string, G4double> m_GlobalParameterMap;
  std::map<std::string, towerposition> m_TowerPostionMap;

  std::set<G4LogicalVolume *> m_AbsorberLogicalVolSet;
  std::set<G4LogicalVolume *> m_ScintiLogicalVolSet;
  std::set<G4LogicalVolume *> m_SupportLogicalVolSet;
};

#endif
