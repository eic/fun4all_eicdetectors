#ifndef G4DETECTORS_PHG4TRDDETECTOR_H
#define G4DETECTORS_PHG4TRDDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <cmath>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class PHG4TRDDetector : public PHG4Detector
{
public:
 PHG4TRDDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam = "TRD", const int lyr =0 );

 //! destructor
  ~PHG4TRDDetector() override
  {
  }

  //! construct TRD
  void ConstructMe(G4LogicalVolume *world) override;
  
 int IsInTRD(const G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
  

 private:
  PHParameters *m_Params;

  
  G4VPhysicalVolume *Phys ;
  G4VPhysicalVolume *fPhysicsRadiator;
  G4VPhysicalVolume *TRD_det_Phys;
  G4VPhysicalVolume *MPGD_win_Phys;
  G4VPhysicalVolume *Cathode_Phys;
  G4VPhysicalVolume *Gas_Active;
  G4VPhysicalVolume *GEM_top_Phys;
  G4VPhysicalVolume *GEM_diel_Phys;
  G4VPhysicalVolume *GEM_bottom_Phys;
  G4VPhysicalVolume *MMG_mesh_Phys;
  G4VPhysicalVolume *Res_lay_Phys;
  G4VPhysicalVolume *MMG_strips_Phys;
  G4VPhysicalVolume *PCB_Phys;
  
  int m_Active;
  int m_AbsorberActive;
  
  int m_Layer;
  std::string m_SuperDetector;
  /*
protected:  
  int m_Active;
  int m_AbsorberActive;
  */

};

#endif
