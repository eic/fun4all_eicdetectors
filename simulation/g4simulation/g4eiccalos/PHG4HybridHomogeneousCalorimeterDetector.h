// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4HYBRIDHOMOGENEOUSCALORIMETERDETECTOR_H
#define G4DETECTORS_PHG4HYBRIDHOMOGENEOUSCALORIMETERDETECTOR_H

#include "PHG4CrystalCalorimeterDefs.h"

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <string>

class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4HybridHomogeneousCalorimeterDisplayAction;
class PHG4Subsystem;
class PHParameters;

/**
 * \file ${file_name}
 * \brief Module to build crystal calorimeter in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4HybridHomogeneousCalorimeterDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4HybridHomogeneousCalorimeterDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4HybridHomogeneousCalorimeterDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //! check if volume is in this calorimeter
  virtual int IsInCrystalCalorimeter(G4VPhysicalVolume *) const;

  // ----- accessing member variables: ------------

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  int get_DetectorId() const { return m_DetectorId; }
  void DetectorId(const int i) { m_DetectorId = i; }
  void DoFullLightProp(bool doProp) { m_doLightProp = doProp; }

  // ----- additional accessors used by derived classes: ------------

  PHParameters *GetParams() { return m_Params; }

 protected:  // for variables also used in PHG4ProjCrystalCalorimeterDetector
  PHG4HybridHomogeneousCalorimeterDisplayAction *GetDisplayAction() { return m_DisplayAction; }
  G4Material *GetCarbonFiber();
  virtual int GetCaloType() const { return PHG4CrystalCalorimeterDefs::CaloType::nonprojective; }

 private:  // private stuff
  G4LogicalVolume* ConstructSupportFrame(G4LogicalVolume *envelope);
  G4LogicalVolume *ConstructTower();
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower);//, G4LogicalVolume *support);
  G4Material* GetScintillatorMaterial(float setting);
  G4Material* GetTedlarMaterial();
  G4Material* GetVM2000Material();
  int ParseParametersFromTable();
  void CrystalTable(G4Material *mat);
  void SurfaceTable(G4LogicalVolume *vol);
  void MakeBoundary(G4VPhysicalVolume *vol1, G4VPhysicalVolume *vol2);
  struct towerposition
  {
    G4double x;
    G4double y;
    G4double z;
    int idx_j;
    int idx_k;
  };

  int m_DetectorId = 0;

  std::string m_SuperDetector;

  PHParameters *m_Params = nullptr;

  PHG4HybridHomogeneousCalorimeterDisplayAction *m_DisplayAction = nullptr;

  std::string _towerlogicnameprefix;

  std::map<std::string, G4double> _map_global_parameter;
  std::map<std::string, towerposition> _map_tower;
  std::set<G4VPhysicalVolume *> m_ActiveVolumeSet;
  std::set<G4VPhysicalVolume *> m_PassiveVolumeSet;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActive;
  int m_AbsorberActive;
  bool m_doLightProp;
};

#endif
