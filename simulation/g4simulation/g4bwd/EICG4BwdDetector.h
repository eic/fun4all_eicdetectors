// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4BwdDETECTOR_H
#define EICG4BwdDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class EICG4BwdDetector : public PHG4Detector
{
 public:
  //! constructor
  EICG4BwdDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int layer = 0);

  //! destructor
  virtual ~EICG4BwdDetector() override {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  //@}

  int GetDetId(G4VPhysicalVolume *) const;
  //void SetSteppingAction(EICG4B0SteppingAction *stpact) { m_SteppingAction = stpact; }
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
  PHParameters *getParams();
  void SetTowerMappingFile(const std::string &filename)
  {
    _mapping_tower_file = filename;
  }
 private:
  G4LogicalVolume *ConstructTower();
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower);
  int ParseParametersFromTable();

  struct towerposition
  {
    G4double x;
    G4double y;
    G4double z;
    int idx_j;
    int idx_k;
  };

//  std::map<std::string, G4double> m_GlobalParameterMap;
  std::map<std::string, towerposition> m_TowerPositionMap;

//  EICG4B0SteppingAction *m_SteppingAction;
  PHParameters *m_Params;
  // active volumes
  std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;
  //  std::set<G4LogicalVolume *>   m_LogicalVolumesSet;
  std::map<G4VPhysicalVolume *, int> m_PhysicalVolumesDet;
  std::set<G4LogicalVolume *> m_LogicalVolSet;
  //  std::map<G4LogicalVolume *, int>   m_LogicalVolumesDet;
  int m_Layer;
  std::string m_SuperDetector;
  std::string _mapping_tower_file;
  std::string m_TowerLogicNamePrefix;
 protected:
 void LogicalVolSetInsert(G4LogicalVolume *logvol)
{
	m_LogicalVolSet.insert(logvol);
}
};

#endif  // EICG4BwdDETECTOR_H
