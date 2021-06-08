// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4ZDCDETECTOR_H
#define EICG4ZDCDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class EICG4ZDCDetector : public PHG4Detector
{
 public:
  //! constructor
  EICG4ZDCDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~EICG4ZDCDetector() {}

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  //@}
  int GetVolumeInfo(G4VPhysicalVolume *volume);

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

 private:
  PHParameters *m_Params;

  std::set<G4LogicalVolume *> m_ActiveLogicalVolumesSet;
  std::set<G4LogicalVolume *> m_AbsorberLogicalVolumesSet;

  // active volumes
  std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;
  std::map<G4LogicalVolume*, int> m_ActiveLogicalVolumeInfoMap;

  std::string m_SuperDetector;
};

#endif // EICG4ZDCDETECTOR_H
