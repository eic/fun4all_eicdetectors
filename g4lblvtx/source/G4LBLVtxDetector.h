// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef G4DETECTORS_G4LBLVTXDETECTOR_H
#define G4DETECTORS_G4LBLVTXDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>

class G4LBLVtxDisplayAction;

class G4LogicalVolume;
class G4VPhysicalVolume;

class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

/*!
 * \brief G4LBLVtxDetector is a generic detector built from a GDML import
 */
class G4LBLVtxDetector : public PHG4Detector
{
 public:
  G4LBLVtxDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam, PHParameters* parameters);

  virtual ~G4LBLVtxDetector();

  //! construct
  void ConstructMe(G4LogicalVolume* world);

  int IsInDetector(G4VPhysicalVolume*) const;

  void SuperDetector(const std::string& name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  void Print(const std::string& what = "ALL") const;

 private:
  void SetActiveVolumes(G4VPhysicalVolume* physvol);

  G4LBLVtxDisplayAction* m_DisplayAction;

  std::string m_GDMPath;
  std::string m_TopVolName;
  std::set<G4VPhysicalVolume*> m_ActivePhysVolumeMap;
  std::set<G4VPhysicalVolume*> m_PassivePhysVolumeMap;
  std::set<std::string> m_ActiveVolName;
  double m_placeX;
  double m_placeY;
  double m_placeZ;

  double m_rotationX;
  double m_rotationY;
  double m_rotationZ;

  int m_Active;
  int m_AbsorberActive;

  std::string m_SuperDetector;
};

#endif /* G4DETECTORS_G4LBLVTXDETECTOR_H */
