// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4RPDETECTOR_H
#define EICG4RPDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class EICG4RPDetector : public PHG4Detector
{
 public:
  //! constructor
  EICG4RPDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int layer = 0);

  //! destructor
  virtual ~EICG4RPDetector() override {}

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  //@}

  int GetDetId(G4VPhysicalVolume *) const;

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
  void SetParametersFromFile();

  PHParameters *getParams();

 private:
  PHParameters *m_Params;
  
  // active volumes (i.e. G4_Si)
  std::map<G4VPhysicalVolume *, int> m_ActivePhysicalVolumesMap;
  // passive volumes (i.e. G4_Cu)
  std::set<G4VPhysicalVolume *> m_PassivePhysicalVolumesSet;
  
  int m_Layer;
  std::string m_SuperDetector;
};

#endif  // EICG4RPDETECTOR_H
