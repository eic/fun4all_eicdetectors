// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERSTRIPDETECTOR_H
#define G4DETECTORS_PHG4CYLINDERSTRIPDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class PHG4CylinderStripDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4CylinderStripDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int layer = 0);

  //! destructor
  virtual ~PHG4CylinderStripDetector()
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume *world);

  //bool IsInDetector(G4VPhysicalVolume *) const;
  bool IsInTileC(G4VPhysicalVolume *) const;
  bool IsInTileZ(G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
  void BuildMaterials();

 private:
  PHParameters *m_Params;

  G4VPhysicalVolume *m_CylinderPhysicalVolume;
  std::set<G4VPhysicalVolume*> m_CylinderCPhysicalVolume;
  std::set<G4VPhysicalVolume*> m_CylinderZPhysicalVolume;
  
  //std::set<G4VPhysicalVolume*> m_PhysicalVolumesSet;

  int m_Layer;
  std::string m_SuperDetector;
  G4VSolid* GetHollowBar();
  double thickness;
  double barwidth;
};

#endif
