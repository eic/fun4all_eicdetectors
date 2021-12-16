// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4LFHCALDETECTOR_H
#define G4DETECTORS_PHG4LFHCALDETECTOR_H

#include <g4main/PHG4Detector.h>
#include <Geant4/G4Material.hh>

#include <Geant4/G4Types.hh>  // for G4double

#include <map>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4LFHcalDisplayAction;
class PHG4Subsystem;
class PHParameters;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4LFHcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4LFHcalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4LFHcalDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //!@name volume accessors
  int IsInLFHcal(G4VPhysicalVolume *) const;

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  PHParameters *getParamsDet() const { return m_Params; }

  int get_Layer() const { return m_Layer; }

  void DoFullLightProp(bool doProp) { m_doLightProp = doProp; }

 private:
  G4LogicalVolume *ConstructTower();
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower);
  int ParseParametersFromTable();
  G4Material *GetScintillatorMaterial();
  G4Material *GetCoatingMaterial();
  G4Material *GetWLSFiberMaterial();
  void SurfaceTable(G4LogicalVolume *vol);
  void MakeBoundary(G4VPhysicalVolume *crystal, G4VPhysicalVolume *opdet, bool isFiber);
  void MakeBoundary_Fiber_Scint(G4VPhysicalVolume *crystal, G4VPhysicalVolume *opdet);
  struct towerposition
  {
    G4double x;
    G4double y;
    G4double z;
    int idx_j;
    int idx_k;
  };

  PHG4LFHcalDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;

  int m_ActiveFlag = 1;
  int m_AbsorberActiveFlag = 0;
  int m_Layer = 0;

  std::string m_TowerLogicNamePrefix;
  std::string m_SuperDetector;

  std::map<std::string, G4double> m_GlobalParameterMap;
  std::map<std::string, towerposition> m_TowerPostionMap;

  std::set<G4LogicalVolume *> m_AbsorberLogicalVolSet;
  std::set<G4LogicalVolume *> m_ScintiLogicalVolSet;

  bool m_doLightProp = false;
};

#endif
