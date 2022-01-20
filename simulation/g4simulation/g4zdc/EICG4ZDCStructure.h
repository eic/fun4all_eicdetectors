#ifndef EICG4ZDCSTRUCTURE_H
#define EICG4ZDCSTRUCTURE_H

#include <set>
#include <map>
#include <Geant4/globals.hh>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VisAttributes;
class G4Material;

class EICG4ZDCStructure {
  public: 
  EICG4ZDCStructure();
  ~EICG4ZDCStructure();

  double ConstructCrystalTowers(double x0, double y0, double z0,
				double x1, double y1, double z1,
				G4VPhysicalVolume *mPhy);

  double ConstructEMLayers(double x0, double y0, double z0, 
			 double x1, double y1, double z1,
			 G4VPhysicalVolume *mPhy);
  double ConstructHCSiliconLayers(double x0, double y0, double z0, 
			      double x1, double y1, double z1,
			      G4VPhysicalVolume *mPhy);
  double ConstructHCSciLayers(double x0, double y0, double z0, 
			      double x1, double y1, double z1,
			      G4VPhysicalVolume *mPhy);
  void ProvideLogicalVolumesSets(std::set<G4LogicalVolume *> &ActiveLogicalVolumesSet,
				 std::set<G4LogicalVolume *> &AbsorberLogicalVolumesSet);
  void ProvideLogicalVolumeInfoMap(std::map<G4LogicalVolume *, int> &ActiveLogicalVolumeInfoMap,
				   std::map<G4LogicalVolume *, int> &AbsorberLogicalVolumeInfoMap);

  void Print();
  void PrintTowerMap(const std::string &d);


private:

  void SetColors();
  void Materials();
  
  int fLayer;
  
  G4Material* fmat_World;
  G4Material* fmat_W;
  G4Material* fmat_PET;
  G4Material* fmat_Sci;
  G4Material* fmat_Si;
  G4Material* fmat_Pb;
  G4Material* fmat_Cu;
  G4Material* fmat_Fe;
  G4Material* fmat_Crystal;

  G4VisAttributes* fvisCrystal;
  G4VisAttributes* fvisPIX;
  G4VisAttributes* fvisPAD;
  G4VisAttributes* fvisDM;
  G4VisAttributes* fvisW;
  G4VisAttributes* fvisPb;
  G4VisAttributes* fvisSci;

  std::set<G4LogicalVolume *> m_ActiveLogicalVolumesSet;
  std::set<G4LogicalVolume *> m_AbsorberLogicalVolumesSet;
  std::map<G4LogicalVolume*, int> m_ActiveLogicalVolumeInfoMap;
  std::map<G4LogicalVolume*, int> m_AbsorberLogicalVolumeInfoMap;

  double _z_Crystal[2];
  double _z_EMLayers[2];
  double _z_HCSilicon[2];
  double _z_HCSci[2];

};  

#endif
