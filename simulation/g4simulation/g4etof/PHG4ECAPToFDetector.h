#ifndef G4DETECTORS_PHG4ECAPToFDETECTOR_H
#define G4DETECTORS_PHG4ECAPToFDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <cmath>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class PHG4ECAPToFDetector : public PHG4Detector
{
 public:
  PHG4ECAPToFDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam = "ETOF", const int lyr = 0);

  //!destructor
  ~PHG4ECAPToFDetector() override
  {
  }

  //Construct mRPC ToF
  void ConstructMe(G4LogicalVolume *world) override;

  int IsInToF(const G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }

 private:
  PHParameters *m_Params;

  G4VPhysicalVolume *Phys;
  G4VPhysicalVolume *fhc_phys;
  G4VPhysicalVolume *fpcb_phys;
  G4VPhysicalVolume *fpcbcu_phys;
  G4VPhysicalVolume *fmylar_phys;
  G4VPhysicalVolume *fcarbon_phys;
  G4VPhysicalVolume *fglass_phys[6];
  G4VPhysicalVolume *fgas_phys[7];
  G4VPhysicalVolume *mcarbon_phys;
  G4VPhysicalVolume *mmylar_phys;
  G4VPhysicalVolume *mpcbcu_phys;
  G4VPhysicalVolume *mpcb_phys;
  G4VPhysicalVolume *mpcbcu2_phys;
  G4VPhysicalVolume *mmylar2_phys;
  G4VPhysicalVolume *mcarbon2_phys;
  G4VPhysicalVolume *bglass_phys[6];
  G4VPhysicalVolume *bgas_phys[7];
  G4VPhysicalVolume *bhc_phys;
  G4VPhysicalVolume *bpcb_phys;
  G4VPhysicalVolume *bpcbcu_phys;
  G4VPhysicalVolume *bmylar_phys;
  G4VPhysicalVolume *bcarbon_phys;

  int m_Active;
  int m_Layer;

  std::string m_SuperDetector;
};
#endif
