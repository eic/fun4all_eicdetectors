// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BEASTMAGNETDETECTOR_H
#define BEASTMAGNETDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>  // for string

class BeastMagnetDisplayAction;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class BeastMagnetDetector : public PHG4Detector
{
 public:
  //! constructor
  BeastMagnetDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~BeastMagnetDetector() {}

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

 private:
  void InsertVolumes(G4VPhysicalVolume* physvol);
BeastMagnetDisplayAction *m_DisplayAction;
  PHParameters *m_Params;

  std::string m_GDMPath;
  std::string m_TopVolName;

  // active volumes
  std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;

  std::string m_SuperDetector;
  int m_AbsorberActive;
};

#endif // BEASTMAGNETDETECTOR_H
