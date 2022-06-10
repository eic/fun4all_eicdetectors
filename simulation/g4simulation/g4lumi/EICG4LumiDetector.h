// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4LUMIDETECTOR_H
#define EICG4LUMIDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <phool/recoConsts.h> //For rc WorldMaterial
#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <Geant4/G4Color.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4TwoVector.hh>      // for G4ThreeVector
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <TSystem.h>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4UniformMagField.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4NistManager.hh>

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <utility>
#include <map>  //?
#include <set> //?
#include <string>  // for string

#include <Geant4/G4CSGSolid.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4NistManager.hh>


class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class EICG4LumiDetector : public PHG4Detector
{
 public:
  //! constructor
  EICG4LumiDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int layer =0);

  //! destructor
  virtual ~EICG4LumiDetector() override {}

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  int IsInVirtualDetector(G4VPhysicalVolume *) const;

  //@}

  int GetDetId(G4VPhysicalVolume *) const;

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }

  void SetParametersFromFile();

  void AddLumiWindow(G4ThreeVector size, G4ThreeVector pos, double angle, std::string material, G4LogicalVolume *logicWorld);
  G4ThreeVector AddLumiExitWindow(G4ThreeVector size, G4ThreeVector pos, double factor, double angle, std::string material, G4LogicalVolume *logicWorld);
  void AddLumiMag_OuterVessel( G4ThreeVector size, G4ThreeVector pos, std::string material, G4LogicalVolume *logicWorld);
  void AddLumiMag_MagCore(G4ThreeVector size, G4ThreeVector pos, double Bx, std::string material, G4LogicalVolume *logicWorld);
  void AddTriangularTrapezoid(G4ThreeVector size, G4ThreeVector pos, double angle, std::string material, G4LogicalVolume *logicWorld);
  void AddCuboid(G4ThreeVector Wsize, G4ThreeVector Wpos, G4ThreeVector Msize, G4ThreeVector Mpos, double angle, std::string material, G4LogicalVolume *logicWorld);
  G4LogicalVolume* AddRecCone(G4ThreeVector Wsize, G4ThreeVector Wpos, G4ThreeVector Msize, G4ThreeVector Mpos, G4ThreeVector Tr2size, G4ThreeVector Tr2pos, double angle, std::string material, G4LogicalVolume *logicWorld);
  void AddExitWindowForV2(G4ThreeVector size, G4ThreeVector pos_daug, double angle, std::string material, G4LogicalVolume *logicRecCone);
  void AddLumiTracker( std::string name, G4ThreeVector size, G4ThreeVector pos, G4LogicalVolume *logicRecCone);
  void AddCAL( std::string name, G4ThreeVector size, G4ThreeVector pos, int total_tower, G4LogicalVolume *logicWorld );
  G4LogicalVolume* MakeTower(G4double calorSizeXY, G4double calorEMZ);
  void AddExitWindowForV3(G4ThreeVector Wsize, G4ThreeVector Wpos, G4ThreeVector Tr2size, G4ThreeVector Tr2pos, std::string material, G4LogicalVolume *logicWorld);

  PHParameters *getParams();

 private:

  PHParameters *m_Params;
   
    // active volumes (e.g. G4_Si)
    std::set<G4VPhysicalVolume *> m_ActivePhysicalVolumesSet;
    // virtual volumes (e.g. G4_Galactic)
    std::map<G4VPhysicalVolume *, int> m_VirtualPhysicalVolumesMap;
    // passive volumes
    std::set<G4VPhysicalVolume *> m_PassivePhysicalVolumesSet; 
    
  int m_Layer;
  std::string m_SuperDetector;

  std::string m_Name;

};

#endif // EICG4LUMIDETECTOR_H
