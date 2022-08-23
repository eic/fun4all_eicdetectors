#include "PHG4FSTDetector.h"
#include "PHG4FSTDisplayAction.h"
#include "PHG4FSTSteppingAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phparameter/PHParameters.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Trd.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>              // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>               // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4PVParameterised.hh>
#include <Geant4/G4PVReplica.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4LogicalSkinSurface.hh>
#include <Geant4/G4LogicalBorderSurface.hh>

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________________
PHG4FSTDetector::PHG4FSTDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters *parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4FSTDisplayAction*>(subsys->GetDisplayAction()))
  , m_SteppingAction(0)
  , _place_in_x(0.0 * mm)
  , _place_in_y(0.0 * mm)
  , _place_in_z(0.0 * mm)
  , _center_offset_x(0.0 * mm)
  , _center_offset_y(0.0 * mm)
  , _quadratic_detector(0)
  , _rot_in_x(0.0)
  , _rot_in_y(0.0)
  , _rot_in_z(0.0)
  , _rMin1(50 * mm)
  , _rMax1(2620 * mm)
  , _rMin2(50 * mm)
  , _rMax2(3369 * mm)
  , _dZ(1000 * mm)
  , _sPhi(0)
  , _dPhi(2 * M_PI)
  , _tower_type(0)
  , _tower_readout(0.5 * mm)
  , _tower_dx(100 * mm)
  , _tower_dy(100 * mm)
  , _tower_dz(1000.0 * mm)
  , _scintFiber_diam(1.0 * mm)
  , _cerenkovFiber_diam(1.0 * mm)
  , _cerenkovFiber_material(0)
  , _tower_makeNotched(0)
  , _absorber_Material(0)
  , _materialScintillator("G4_POLYSTYRENE")
  , _materialAbsorber("G4_Fe")
  , _active(1)
  , _absorberactive(0)
  , _layer(0)
  , _blackhole(0)
  , _towerlogicnameprefix("hFSTTower")
  , _superdetector("NONE")
  , _mapping_tower_file("")
  , m_Params(parameters)
{
}
//_______________________________________________________________________
int PHG4FSTDetector::IsInActiveSensorFST(G4VPhysicalVolume* volume,const std::string supdet) const
{
  if (volume->GetName().find(supdet +"_SiliconSensor") != string::npos)
  {
    if (_active)
      return 1;
    else
      return 0;
  }
  //only record energy in actual absorber- drop energy lost in air gaps inside FST envelope
  else if (volume->GetName().find("absorber") != string::npos)
  {
    if (_absorberactive)
      return -1;
    else
      return 0;
  }
  else if (volume->GetName().find("envelope") != string::npos)
  {
    return 0;
  }
  return 0;
}

//_______________________________________________________________________
void PHG4FSTDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4FSTDetector: Begin Construction" << endl;
  }

    ConstructSTDisk(logicWorld);
  return;
}

void PHG4FSTDetector::ConstructSTDisk(G4LogicalVolume* mother){


  G4double place_z = m_Params->get_double_param("z_position") *cm;
  G4double min_radius = m_Params->get_double_param("r_min") * cm;
  G4double max_radius = m_Params->get_double_param("r_max") *cm;
  G4double silicon_thickness = 35 *um;
  G4double offset_cutout = m_Params->get_double_param("offset_cutout")*cm;

  string  layer_name[] = {"SiliconSensor", "Metalconnection", "HDI", "Cooling", "Support", "Support_Gap", "Support2"};
  G4Material* materialLayer[] = {GetDetectorMaterial("G4_Si"), GetDetectorMaterial("G4_Al"), GetDetectorMaterial("G4_KAPTON"), GetDetectorMaterial("G4_WATER"), GetDetectorMaterial("G4_GRAPHITE"), GetDetectorMaterial("G4_AIR"), GetDetectorMaterial("G4_GRAPHITE")};
  G4double thicknessLayer[] = {silicon_thickness, 15 * um, 20 * um, 100 * um, 50 * um, 1 * cm, 50 * um};

  G4double totalThickness = 0;
  for(int ilay=0; ilay<7; ilay++){
    totalThickness += thicknessLayer[ilay];
  }
  G4double currentPosition = 0;
  for(int ilay=0; ilay<7; ilay++){
    G4VSolid *sol_layer = new G4Tubs("sol_tmp_" + layer_name[ilay],
                                                0,
                                                max_radius,
                                                thicknessLayer[ilay] / 2,
                                                0.,2*M_PI);
    G4VSolid *sol_layer_cutout = new G4Tubs("sol_cutout_" + layer_name[ilay],
                                                0,
                                                min_radius,
                                                thicknessLayer[ilay],
                                                0.,2*M_PI);
    sol_layer = new G4SubtractionSolid("sol_" + layer_name[ilay], sol_layer, sol_layer_cutout, 0, G4ThreeVector(offset_cutout,0,0));

    G4LogicalVolume *Log_layer = new G4LogicalVolume(sol_layer, materialLayer[ilay], "log_" + layer_name[ilay]);
    m_DisplayAction->AddVolume(Log_layer, layer_name[ilay]);

    currentPosition += thicknessLayer[ilay] / 2;
    if(place_z<0){
      new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, place_z+totalThickness/2-currentPosition),
                            Log_layer,
                            _superdetector+"_"+layer_name[ilay]+std::to_string(ilay),
                            mother,
                            0, 0, OverlapCheck());
    } else {
      new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, place_z-totalThickness/2+currentPosition),
                            Log_layer,
                            _superdetector+"_"+layer_name[ilay]+std::to_string(ilay),
                            mother,
                            0, 0, OverlapCheck());
    }
    currentPosition += thicknessLayer[ilay] / 2;
  }

  return;
}

G4Material* PHG4FSTDetector::MakeCarbonFoamMaterial_Longeron(){
  G4Material* carbon_foam = GetDetectorMaterial("C_FOAM_FST_LONGERON", false);  // false suppresses warning that material does not exist
  if(!carbon_foam){
    G4double density;
    G4int ncomponents;
    carbon_foam = new G4Material("C_FOAM_FST_LONGERON", density = 0.07 * g / cm3, ncomponents = 2); // VERY CONSERVATIVE DENSITY
    // carbon_foam = new G4Material("C_FOAM_FST", density = 0.26 * g / cm3, ncomponents = 2); // CONSERVATIVE DENSITY
    // carbon_foam = new G4Material("C_FOAM_FST", density = 0.06 * g / cm3, ncomponents = 2); // LIGHTEST DENSITY
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.97);
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 0.03);
  }
  return carbon_foam;

}

G4Material* PHG4FSTDetector::MakeCarbonFoamMaterial_Wheel(){
  G4Material* carbon_foam = GetDetectorMaterial("C_FOAM_FST_WHEEL", false);  // false suppresses warning that material does not exist
  if(!carbon_foam){
    G4double density;
    G4int ncomponents;
    carbon_foam = new G4Material("C_FOAM_FST_WHEEL", density = 0.20 * g / cm3, ncomponents = 2); // VERY CONSERVATIVE DENSITY
    // carbon_foam = new G4Material("C_FOAM_FST", density = 0.26 * g / cm3, ncomponents = 2); // CONSERVATIVE DENSITY
    // carbon_foam = new G4Material("C_FOAM_FST", density = 0.06 * g / cm3, ncomponents = 2); // LIGHTEST DENSITY
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.97);
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 0.03);
  }
  return carbon_foam;

}
G4Material* PHG4FSTDetector::MakeCarbonFleece(){
  G4Material* carbon_foam = GetDetectorMaterial("C_FLEECE_FST", false);  // false suppresses warning that material does not exist
  if(!carbon_foam){
    G4double density;
    G4int ncomponents;
    carbon_foam = new G4Material("C_FLEECE_FST", density = 0.40 * g / cm3, ncomponents = 2); // VERY CONSERVATIVE DENSITY
    // carbon_foam = new G4Material("C_FOAM_FST", density = 0.26 * g / cm3, ncomponents = 2); // CONSERVATIVE DENSITY
    // carbon_foam = new G4Material("C_FOAM_FST", density = 0.06 * g / cm3, ncomponents = 2); // LIGHTEST DENSITY
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.97);
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 0.03);
  }
  return carbon_foam;

}

G4Material* PHG4FSTDetector::GetCarbonFiber()
{
  static string matname = "FSTCarbonFiber";
  G4Material* carbonfiber = G4Material::GetMaterial(matname, false);  // false suppresses warning that material does not exist
  if (!carbonfiber)
  {
    G4double density_carbon_fiber = 1.44 * g / cm3;
    carbonfiber = new G4Material(matname, density_carbon_fiber, 1);
    carbonfiber->AddElement(G4Element::GetElement("C"), 1);
  }
  return carbonfiber;
}

G4Material* PHG4FSTDetector::GetKapton()
{
  static string matname = "FSTKapton";
  G4Material* kaptonmat = G4Material::GetMaterial(matname, false);  // false suppresses warning that material does not exist
  if (!kaptonmat)
  {
    G4double density_carbon_fiber = 2.02 * g / cm3; // 28.41cm rad length
    G4int ncomponents;
    kaptonmat = new G4Material(matname, density_carbon_fiber, ncomponents = 4);
    kaptonmat->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 22);
    kaptonmat->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 10);
    kaptonmat->AddElement(G4NistManager::Instance()->FindOrBuildElement("N"), 2);
    kaptonmat->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 5);
  }
  return kaptonmat;
}


G4Material* PHG4FSTDetector::MakeCarbonHoneyCombMaterial()
{
  G4Material* carbonfiber = G4Material::GetMaterial("TTLCarbonHoneyComb", false);  // false suppresses warning that material does not exist
  if (!carbonfiber)
  {
    G4double density_carbon_fiber = 0.03 * g / cm3;
    carbonfiber = new G4Material("TTLCarbonHoneyComb", density_carbon_fiber, 3);
    carbonfiber->AddElement(G4Element::GetElement("O"), 0.074);
    carbonfiber->AddElement(G4Element::GetElement("C"), 0.903);
    // carbonfiber->AddElement(G4Element::GetElement("C"), 0.870);
    carbonfiber->AddElement(G4Element::GetElement("H"), 0.023);
    // carbonfiber->AddElement(G4Element::GetElement("G4_Cl"), 0.033);
  }
  return carbonfiber;
}

G4Material* PHG4FSTDetector::MakeGlue()
{
  G4Material* carbonfiber = G4Material::GetMaterial("FST_GLUE", false);  // false suppresses warning that material does not exist
  if (!carbonfiber)
  {
    G4double density_carbon_fiber = 0.8 * g / cm3;
    carbonfiber = new G4Material("FST_GLUE", density_carbon_fiber, 3);
    carbonfiber->AddElement(G4Element::GetElement("O"), 0.074);
    carbonfiber->AddElement(G4Element::GetElement("C"), 0.903);
    // carbonfiber->AddElement(G4Element::GetElement("C"), 0.870);
    carbonfiber->AddElement(G4Element::GetElement("H"), 0.023);
    // carbonfiber->AddElement(G4Element::GetElement("G4_Cl"), 0.033);
  }
  return carbonfiber;
}

