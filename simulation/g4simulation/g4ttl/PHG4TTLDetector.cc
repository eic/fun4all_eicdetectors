#include "PHG4TTLDetector.h"
#include "PHG4TTLDisplayAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>      // for PHG4Subsystem

#include <phparameter/PHParameters.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PVReplica.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Trd.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4double, G4int
#include <Geant4/G4AssemblyVolume.hh>  // for G4double, G4int
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4RotationMatrix.hh>

#include <algorithm>  // for max
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <map>  // for _Rb_tree_iterator, _Rb_tree_co...
#include <sstream>
#include <utility>  // for pair

class G4VPhysicalVolume;
class PHCompositeNode;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4TTLDetector::PHG4TTLDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , name_base(dnam)
  , m_DisplayAction(dynamic_cast<PHG4TTLDisplayAction *>(subsys->GetDisplayAction()))
  , m_SteppingAction(0)
  , m_Params(parameters)
  , overlapcheck_sector(false)
{
}

//_______________________________________________________________
//_______________________________________________________________
bool PHG4TTLDetector::IsInSectorActive(G4VPhysicalVolume *physvol)
{
  for (map_phy_vol_t::const_iterator it = map_active_phy_vol.begin();
       it != map_active_phy_vol.end(); ++it)
  {
    if (physvol == (*it).second)
    {
      return true;
    }
  }

  return false;
}

//_______________________________________________________________
void PHG4TTLDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  if (m_Params->get_int_param("isForward") == 1)
  {
    BuildForwardTTL(logicWorld);
  }
  else
  {
    if(m_Params->get_int_param("isECCE"))
      BuildBarrelTTL(logicWorld);
    else
      BuildBarrelTTLStaves(logicWorld);
  }
}

void PHG4TTLDetector::BuildBarrelTTLStaves(G4LogicalVolume *logicWorld)
{

  G4double rCenter = m_Params->get_double_param("rMin");  // center location of Al support plate
  // G4double det_height = 2.1 * cm;
  G4double place_z = m_Params->get_double_param("place_z");
  // G4ThreeVector detzvec(0, 0, place_z);
  G4double detlength = m_Params->get_double_param("length");

  G4double carbonsupport_height = 0.6 * cm;
  G4double carbonfoam_length = 3.1 * cm;
  G4double carbonhoneycomb_length = 2.5 * cm;

  G4double sensor_width = 3.2 * cm;
  G4double sensor_length = 5.5 * cm;
  G4double sensor_thickness = 0.3 * mm;
  G4double sensor_overlap = 3.0 * mm;

  G4double ASIC_thickness = 0.3 * mm;
  G4double ASIC_width = 1.0 * cm;

  G4double module_width = 5.6 * cm;
  G4double module_length = detlength ; //67.5 * cm;
  G4double module_height = sensor_thickness>ASIC_thickness ? carbonsupport_height+sensor_thickness : carbonsupport_height+ASIC_thickness;

  // make module mother volume
  G4VSolid *sol_module_box = new G4Box("sol_module_box",
                                          module_width / 2,
                                          module_height / 2,
                                          module_length / 2);
  G4LogicalVolume *log_module_box = new G4LogicalVolume(sol_module_box, GetDetectorMaterial("G4_AIR"), "log_module_box");
  RegisterLogicalVolume(log_module_box);
  m_DisplayAction->AddVolume(log_module_box, "StripBox");

  // create individual module components
  // carbon honeycomb with cooling pipe
  G4VSolid *sol_carbonhoneycomb_module_mother = new G4Box("sol_carbonhoneycomb_module_mother",
                                          carbonhoneycomb_length / 2,
                                          carbonsupport_height / 2,
                                          module_length / 2);
  G4VSolid *sol_carbonhoneycomb_module = new G4Box("sol_carbonhoneycomb_module_tmp",
                                          carbonhoneycomb_length / 2,
                                          carbonsupport_height / 2,
                                          module_length / 2);

  G4double diameter_coolingtube = 5 * mm;
  G4double wallthickness_coolingtube = 1 * mm;
  G4VSolid *sol_cutout_tube = new G4Tubs("sol_cutout_tube",
                                          0.,
                                          1.003*diameter_coolingtube / 2,
                                          (module_length * 2) / 2,
                                          0.,2*M_PI);
  G4VSolid *sol_cooling_tube = new G4Tubs("sol_cooling_tube_tmp",
                                          (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                          diameter_coolingtube / 2,
                                          (module_length - 0.2 * mm) / 2,
                                          0.,2*M_PI);
  sol_carbonhoneycomb_module = new G4SubtractionSolid(G4String("sol_carbonhoneycomb_module"), sol_carbonhoneycomb_module, sol_cutout_tube, 0, G4ThreeVector(0,0,0));

  G4LogicalVolume *Log_cooling_tube = new G4LogicalVolume(sol_cooling_tube, GetDetectorMaterial("G4_Al"), "Log_cooling_tube");
  RegisterLogicalVolume(Log_cooling_tube);
  m_DisplayAction->AddVolume(Log_cooling_tube, "Cooling_tube");

  G4VSolid *sol_water_cooling = new G4Tubs("sol_water_cooling",
                                          0,
                                          0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                          (module_length - 0.2 * mm) / 2,
                                          0.,2*M_PI);
  G4LogicalVolume *Log_water_cooling = new G4LogicalVolume(sol_water_cooling,  //
                                                           GetDetectorMaterial("G4_WATER"), "Log_water_cooling");
  RegisterLogicalVolume(Log_water_cooling);
  m_DisplayAction->AddVolume(Log_water_cooling, "Water_cooling");
  
  
  G4LogicalVolume *log_carbonhoneycomb_module_mother = new G4LogicalVolume(sol_carbonhoneycomb_module_mother, GetDetectorMaterial("G4_AIR"), "log_carbonhoneycomb_module_mother");
  G4LogicalVolume *log_carbonhoneycomb_module = new G4LogicalVolume(sol_carbonhoneycomb_module, MakeCarbonFoamMaterial(), "log_carbonhoneycomb_module");

  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-(carbonhoneycomb_length+carbonfoam_length)/2 + carbonhoneycomb_length/2, module_height/2-carbonsupport_height/2, 0), log_carbonhoneycomb_module_mother,
                                            "physical_carbonhoneycomb_module_mother", log_module_box, false, 0, overlapcheck_sector),false);
  RegisterLogicalVolume(log_carbonhoneycomb_module_mother);
  m_DisplayAction->AddVolume(log_carbonhoneycomb_module_mother, "Invisible");

  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, 0), log_carbonhoneycomb_module,
                                            "physical_carbonhoneycomb_module", log_carbonhoneycomb_module_mother, false, 0, overlapcheck_sector),false);
  RegisterLogicalVolume(log_carbonhoneycomb_module);
  m_DisplayAction->AddVolume(log_carbonhoneycomb_module, "CarbonHoneycomb");

  // place cooling pipe and water
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Log_cooling_tube,
                                            "physical_cooling_tube", log_carbonhoneycomb_module_mother, false, 0, overlapcheck_sector),false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Log_water_cooling,
                                            "physical_water_cooling", log_carbonhoneycomb_module_mother, false, 0, overlapcheck_sector),false);

  // carbon foam
  G4VSolid *sol_carbonfoam_module = new G4Box("sol_carbonfoam_module",
                                          carbonfoam_length / 2,
                                          carbonsupport_height / 2,
                                          module_length / 2);
  G4LogicalVolume *log_carbonfoam_module = new G4LogicalVolume(sol_carbonfoam_module, MakeCarbonFoamMaterial(), "log_carbonfoam_module");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(((carbonfoam_length+carbonhoneycomb_length)/2 - carbonfoam_length/2), module_height/2-carbonsupport_height/2, 0), log_carbonfoam_module,
                                            "physical_carbonfoam_module", log_module_box, false, 0, overlapcheck_sector),false);
  RegisterLogicalVolume(log_carbonfoam_module);
  m_DisplayAction->AddVolume(log_carbonfoam_module, "CarbonFoam");

  // sensor 3.2 x 5.5 cm
  G4double sensor_gap = 0.1 * cm;
  int nSensorsMother = (int) (module_length / (sensor_length+sensor_gap));
  std::cout << "nSensorsMother: " << nSensorsMother << std::endl;
  // module long envelope for sensors
  G4VSolid *sol_sensor_envelope = new G4Box("sol_sensor_envelope",
                                          sensor_width / 2,
                                          sensor_thickness / 2,
                                          module_length / 2);
  G4LogicalVolume *log_sensor_envelope = new G4LogicalVolume(sol_sensor_envelope, GetDetectorMaterial("G4_AIR"), "log_sensor_envelope");
  // element for sensor that gets replicated
  G4VSolid *sol_sensor_mother = new G4Box("sol_sensor",
                                          sensor_width / 2,
                                          sensor_thickness / 2,
                                          (sensor_length+sensor_gap) / 2);
  G4LogicalVolume *log_sensor_mother = new G4LogicalVolume(sol_sensor_mother, GetDetectorMaterial("G4_AIR"), "log_sensor_mother");
  // sensor that is placed in replicated volume
  G4VSolid *sol_sensor = new G4Box("sol_sensor",
                                          sensor_width / 2,
                                          sensor_thickness / 2,
                                          sensor_length / 2);
  G4LogicalVolume *log_sensor = new G4LogicalVolume(sol_sensor, GetDetectorMaterial("G4_Si"), "log_sensor");

  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, 0), log_sensor,
                                            "physical_Sensor", log_sensor_mother, false, 0, overlapcheck_sector),true);

  // place sensors along log_module_box
  new G4PVReplica("placed_sensor_replicas", log_sensor_mother, log_sensor_envelope,
                kZAxis, nSensorsMother, sensor_length+sensor_gap, 0);

  RegisterLogicalVolume(log_sensor);
  RegisterLogicalVolume(log_sensor_envelope);
  RegisterLogicalVolume(log_sensor_mother);
  m_DisplayAction->AddVolume(log_sensor, "Sensor");
  m_DisplayAction->AddVolume(log_sensor_envelope, "Black");
  m_DisplayAction->AddVolume(log_sensor_mother, "Black");


  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(module_width/2-sensor_width/2, -module_height/2 + sensor_thickness/2, 0), log_sensor_envelope,
                                            "physical_SensorReplicas", log_module_box, false, 0, overlapcheck_sector),false);

  // frontend ASIC
  G4VSolid *sol_ASIC = new G4Box("sol_ASIC",
                                          ASIC_width / 2,
                                          ASIC_thickness / 2,
                                          module_length / 2);
  G4LogicalVolume *log_ASIC = new G4LogicalVolume(sol_ASIC, GetDetectorMaterial("G4_Si"), "log_ASIC");


  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(module_width/2-sensor_width - 3*mm - ASIC_width/2, -module_height/2 + ASIC_thickness/2, 0), log_ASIC,
                                            "physical_ASIC", log_module_box, false, 0, overlapcheck_sector),false);
  RegisterLogicalVolume(log_ASIC);
  m_DisplayAction->AddVolume(log_ASIC, "ASIC");

  // place log_module_box at radial position rCenter tilted by 18 degrees
  // make loop to place module at multiple places in azimuth at radius rCenter
  // int nAzimuth = m_Params->get_int_param("nAzimuth");
  // int nAzimuth = 1 + sensor_overlap - sensor_overlap;
  int nAzimuth = 2 * M_PI * rCenter / (sensor_width * cos(18 * deg) - sensor_overlap);
  std::cout << "nAzimuth = " << nAzimuth << std::endl;
  for(int i = 0; i < nAzimuth; i++)
  {
    G4RotationMatrix *rot_module = new G4RotationMatrix();
    rot_module->rotateZ(i * 360.0 * deg / nAzimuth + 18.0 * deg);
    G4double phi = i * 2 * M_PI / nAzimuth;//m_Params->get_int_param("n_azimuthal")
    G4ThreeVector detzvec(rCenter * sin(phi), rCenter * cos(phi), place_z);
    RegisterPhysicalVolume(new G4PVPlacement(rot_module, detzvec, log_module_box,
                                              "physical_module_"+std::to_string(i), logicWorld, false, i, overlapcheck_sector),false);
  }

}

void PHG4TTLDetector::BuildBarrelTTL(G4LogicalVolume *logicWorld)
{
  // MATERIALS
  G4double density, a;
  G4int ncomponents, z;
  G4int natoms;
  G4String symbol;

  m_SteppingAction->SetNPhiModules(12);
  m_SteppingAction->SetIsForwardTTL(false);

  G4NistManager* man = G4NistManager::Instance();
  G4Element *elH = new G4Element("Hydrogen", symbol = "H", 1., 1.01 * g / mole);
  G4Element *elC = new G4Element("Carbon", symbol = "C", 6., 12.01 * g / mole);
  G4Element *elN = new G4Element("Nitrogen", symbol = "N", 7., 14.01 * g / mole);
  G4Element *elO = new G4Element("Oxygen", symbol = "O", 8., 16.00 * g / mole);
  G4Material *mat_Epoxy = GetDetectorMaterial("EpoxyTTL",false);
  if (!mat_Epoxy)
  {
    mat_Epoxy = new G4Material("EpoxyTTL", density = 1.16 * g / cm3, natoms = 4);
    mat_Epoxy->AddElement(elH, 32);  // Hydrogen
    mat_Epoxy->AddElement(elN, 2);   // Nitrogen
    mat_Epoxy->AddElement(elO, 4);   // Oxygen
    mat_Epoxy->AddElement(elC, 15);  // Carbon
    // G4Material *mat_Epoxy = GetDetectorMaterial("Epoxy");
  }
  G4Material *mat_ALN = GetDetectorMaterial("AluminiumNitrate",false);
  if (!mat_ALN)
  {
    mat_ALN = new G4Material("AluminiumNitrate", density = 3.255 * g / cm3, ncomponents = 2);
    // G4Material *mat_ALN = new G4Material("AluminiumNitrate", density = 3.255 * g / cm3, ncomponents = 2);
    mat_ALN->AddElement(G4Element::GetElement("Al"), 1);
    mat_ALN->AddElement(G4Element::GetElement("N"), 1);
  }
  G4Material *mat_Solder_Tin = GetDetectorMaterial("Tin",false);
  if (!mat_Solder_Tin)
  {
    mat_Solder_Tin = new G4Material("Tin", z = 50., a = 118.7 * g / mole, density = 7.310 * g / cm3);
  }
  G4Material *Air = GetDetectorMaterial("G4_AIR");

  // Carbon fiber 190 width + 65 microns scotch , overall 255 microns
  G4Material *CarbonFiber = new G4Material("CarbonFiber",  density =  1.750*g/cm3, ncomponents=2);
  CarbonFiber->AddMaterial(man->FindOrBuildMaterial("G4_C"), 74.5*perCent);  // Carbon
  CarbonFiber->AddMaterial(mat_Epoxy,                           25.50*perCent);  // Epoxy (scotch)

  // positions
  G4double rCenter = m_Params->get_double_param("rMin");  // center location of Al support plate
  G4double det_height = 2.1 * cm;
  G4double place_z = m_Params->get_double_param("place_z");
  G4ThreeVector detzvec(0, 0, place_z);
  G4double detlength = m_Params->get_double_param("length");

  //Create the envelope = 'world volume' for the calorimeter
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  // Single module with length based on readout (contains 14 LGADs [counting across both sides] in x-direction and 6 in z-direction)
  G4double baseplate_length = 43.1 * mm;
  G4double baseplate_width = 56.5 * mm / 2;
  G4double segmentlength = 6 * baseplate_length;  //(detlength - 10 * cm) / 6;//m_Params->get_double_param("length");

  G4VSolid *sol_module_envelope = new G4Trd("sol_module_envelope",
                                            sin(M_PI / 12.) * (rCenter-det_height/2), sin(M_PI / 12.) * (rCenter + det_height/2),
                                            segmentlength / 2, segmentlength / 2,
                                            (det_height*1.001) / 2);

  G4LogicalVolume *log_module_envelope = new G4LogicalVolume(sol_module_envelope, Air, "log_module_envelope");

  G4double diameter_coolingtube = 5 * mm;
  G4double cooling_plate_height = 1 * mm;
  if(m_Params->get_double_param("cooling_plate_height")>0) cooling_plate_height = m_Params->get_double_param("cooling_plate_height");

  G4VSolid *sol_cooling_plate_top = new G4Box("sol_cooling_plate_top",
                                          sin(M_PI / 12.) * (rCenter + diameter_coolingtube / 2 ),
                                          segmentlength / 2,
                                          cooling_plate_height / 2);
  G4VSolid *sol_cooling_plate_bottom = new G4Box("sol_cooling_plate_top",
                                          sin(M_PI / 12.) * (rCenter - diameter_coolingtube / 2 - cooling_plate_height),
                                          segmentlength / 2,
                                          cooling_plate_height / 2);
  // std::cout << "top plate: " << sin(M_PI / 12.) * (rCenter + diameter_coolingtube / 2 + cooling_plate_height / 2) << std::endl;
  // std::cout << "bottom plate: " << sin(M_PI / 12.) * (rCenter - diameter_coolingtube / 2 - cooling_plate_height / 2) << std::endl;
  G4LogicalVolume *log_cooling_plate_top = new G4LogicalVolume(sol_cooling_plate_top, GetDetectorMaterial("G4_Al"), "log_cooling_plate_barrel_top");
  G4LogicalVolume *log_cooling_plate_bottom = new G4LogicalVolume(sol_cooling_plate_bottom, GetDetectorMaterial("G4_Al"), "log_cooling_plate_barrel_bottom");

  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, diameter_coolingtube/2+cooling_plate_height/2), log_cooling_plate_top,
                                            "physical_cooling_plate_top", log_module_envelope, false, 0, overlapcheck_sector),false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, -(diameter_coolingtube/2+cooling_plate_height/2)), log_cooling_plate_bottom,
                                            "physical_cooling_plate_bottom", log_module_envelope, false, 0, overlapcheck_sector),false);
  RegisterLogicalVolume(log_cooling_plate_top);
  RegisterLogicalVolume(log_cooling_plate_bottom);
  m_DisplayAction->AddVolume(log_cooling_plate_top, "CoolingPlate");
  m_DisplayAction->AddVolume(log_cooling_plate_bottom, "CoolingPlate");

  G4double wallthickness_coolingtube = 1 * mm;
  G4VSolid *sol_cutout_tube = new G4Box("sol_cutout_tube",
                                          (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                          (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                          (segmentlength * 1.1) / 2);
  G4VSolid *sol_cooling_tube = new G4Box("sol_cooling_tube_tmp",
                                          diameter_coolingtube / 2,
                                          diameter_coolingtube / 2,
                                          (segmentlength - 0.2 * mm) / 2);
  sol_cooling_tube = new G4SubtractionSolid(G4String("sol_cooling_tube"), sol_cooling_tube, sol_cutout_tube, 0, G4ThreeVector(0,0,0));

  G4LogicalVolume *Log_cooling_tube = new G4LogicalVolume(sol_cooling_tube, GetDetectorMaterial("G4_Al"), "Log_cooling_tube");
  RegisterLogicalVolume(Log_cooling_tube);
  m_DisplayAction->AddVolume(Log_cooling_tube, "Cooling_tube");

  G4VSolid *sol_water_cooling = new G4Box("sol_water_cooling",
                                          0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                          0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                          (segmentlength - 0.2 * mm) / 2);
  G4LogicalVolume *Log_water_cooling = new G4LogicalVolume(sol_water_cooling,  //
                                                           GetDetectorMaterial("G4_WATER"), "Log_water_cooling");
  RegisterLogicalVolume(Log_water_cooling);
  m_DisplayAction->AddVolume(Log_water_cooling, "Water_cooling");

  G4VSolid *sol_internal_support_center = new G4Box("sol_internal_support_center",
                                          (1.5 * baseplate_width - diameter_coolingtube)/2,
                                          diameter_coolingtube / 2,
                                          (2* mm) / 2);

  G4LogicalVolume *Log_internal_support_center = new G4LogicalVolume(sol_internal_support_center, CarbonFiber, "Log_internal_support_center");
  RegisterLogicalVolume(Log_internal_support_center);
  m_DisplayAction->AddVolume(Log_internal_support_center, "Carbon_Support");
  G4VSolid *sol_internal_support_edge = new G4Box("sol_internal_support_edge",
                                          (1.5 * baseplate_width - diameter_coolingtube - baseplate_width / 3)/2,
                                          diameter_coolingtube / 2,
                                          (2* mm) / 2);

  G4LogicalVolume *Log_internal_support_edge = new G4LogicalVolume(sol_internal_support_edge, CarbonFiber, "Log_internal_support_edge");
  RegisterLogicalVolume(Log_internal_support_edge);
  m_DisplayAction->AddVolume(Log_internal_support_edge, "Carbon_Support");


  G4RotationMatrix *rotcooling = new G4RotationMatrix();
  rotcooling->rotateX(M_PI / 2);
  G4double leftedgeCU = sin(M_PI / 12.) * (rCenter + det_height / 2 + cooling_plate_height / 2);
  int maxicup = 9;
  for (int icup = 0; icup < maxicup; icup++)
  {
    G4double edgeshift = 0;
    if (icup == 0) edgeshift = baseplate_width / 3;
    if (icup == (maxicup - 1)) edgeshift = -baseplate_width / 3;
    RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift, 0, 0), Log_cooling_tube,
                                        "cooling_tube_Physical_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),   false);
    RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift, 0, 0), Log_water_cooling,
                                        "water_cooling_Physical_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),   false);
    if(icup!=0 && icup<(maxicup-2)){
      RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift + 0.75*baseplate_width, 0, 0), Log_internal_support_center,
                                        "internal_support_center_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),   false);
      RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift + 0.75*baseplate_width, 2*baseplate_length, 0), Log_internal_support_center,
                                        "internal_support_center_f_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),   false);
      RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift + 0.75*baseplate_width, -2*baseplate_length, 0), Log_internal_support_center,
                                        "internal_support_center_b_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),   false);
    } else if(icup==0 || icup==(maxicup-2)){
      RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift + 0.75*baseplate_width - baseplate_width / 6, 0, 0), Log_internal_support_edge,
                                        "internal_support_edge_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),   false);
      RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift + 0.75*baseplate_width - baseplate_width / 6, 2*baseplate_length, 0), Log_internal_support_edge,
                                        "internal_support_edge_f_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),   false);
      RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift + 0.75*baseplate_width - baseplate_width / 6, -2*baseplate_length, 0), Log_internal_support_edge,
                                        "internal_support_edge_b_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),   false);
    }
  }


  // Sensor Module:
  G4double sensor_width = 21.2 * mm;
  G4double sensor_length = 42.0 * mm;
  G4double baseSH_width_top = baseplate_width / 2 - (0.1543*cm/2);  //-0.15*mm;
  G4double baseSH_width = baseplate_width / 2 - (0.232*cm/2);  //-0.15*mm;

  const int nLayers = 8;
  std::string strLayerName[nLayers] = {
      "ThermalPad",
      "ALN",
      "LairdFilm",
      "ROC",
      "Solder",
      "Sensor",
      "Epoxy",
      "AIN"};
  G4Material *materialLayer[nLayers] = {GetDetectorMaterial("G4_GRAPHITE"),
                                        mat_ALN,
                                        GetDetectorMaterial("G4_GRAPHITE"),
                                        GetDetectorMaterial("G4_PLEXIGLASS"),
                                        mat_Solder_Tin,
                                        GetDetectorMaterial("G4_Si"),
                                        mat_Epoxy,
                                        mat_ALN};
  G4double thicknessLayer[nLayers] = {0.25 * mm, 0.79 * mm, 0.08 * mm, 0.25 * mm, 0.03 * mm, 0.3 * mm, 0.08 * mm, 0.51 * mm};
  G4double widthLayer[nLayers] = {
      baseplate_width - 0.3 * mm,
      baseplate_width,
      sensor_width + 1 * mm,
      sensor_width + 1 * mm,
      sensor_width - 0.2 * mm,
      sensor_width,
      sensor_width,
      baseplate_width - 4 * mm};
  G4double offsetLayer[nLayers] = {
      0.3 * mm / 2,
      0,
      (baseplate_width - widthLayer[2]) / 2 - 0.2 * mm,
      (baseplate_width - widthLayer[3]) / 2 - 0.2 * mm,
      (baseplate_width - widthLayer[4]) / 2 - 0.1 * mm,
      (baseplate_width - widthLayer[5] - 0.1 * mm) / 2,
      (baseplate_width - widthLayer[6] - 0.1 * mm) / 2,
      4 * mm / 2 - 0.2 * mm};
  G4double lengthLayer[nLayers] = {
      baseplate_length - 0.2 * mm,
      baseplate_length,
      sensor_length + 0.2 * mm,
      sensor_length + 0.2 * mm,
      sensor_length - 0.2 * mm,
      sensor_length,
      sensor_length,
      baseplate_length - 0.2 * mm};
  bool layerActive[nLayers] = {
      false, false, false, false, false, true, false, false};

  G4double thicknessDet = 0;
  for (int ilay = 0; ilay < nLayers; ilay++)
  {
    thicknessDet += thicknessLayer[ilay];
  }

  // Sensor Ladder (6 Sensors)
  G4VSolid *sol_sensor_ladder = new G4Box("sol_sensor_ladder",
                                          baseplate_width / 2,
                                          segmentlength / 2,
                                          thicknessDet / 2);
  G4LogicalVolume *log_sensor_ladder = new G4LogicalVolume(sol_sensor_ladder, Air, "log_sensor_ladder");
  m_DisplayAction->AddVolume(log_sensor_ladder, "SensorLadder");

  G4VSolid *sol_sensor_stack = new G4Box("sol_sensor_stack",
                                          baseplate_width / 2,
                                          baseplate_length / 2,
                                          thicknessDet / 2);
  G4LogicalVolume *log_sensor_stack = new G4LogicalVolume(sol_sensor_stack, Air, "log_sensor_stack");
  m_DisplayAction->AddVolume(log_sensor_stack, "SensorStack");

  double z_start = -thicknessDet / 2;
  for (int ilay = 0; ilay < nLayers; ilay++)
  {
    const std::string layer_name = "sensor_stack_" + strLayerName[ilay];
    const std::string layer_name_Solid = "sol_" + layer_name;

    G4VSolid *sol_Module_Layer_Raw = new G4Box(layer_name_Solid + "_Raw",
                                                widthLayer[ilay] / 2,
                                                lengthLayer[ilay] / 2,
                                                thicknessLayer[ilay] / 2);

    G4LogicalVolume *Log_Layer = new G4LogicalVolume(sol_Module_Layer_Raw,  //
                                                     materialLayer[ilay], layer_name + "_Log");
    RegisterLogicalVolume(Log_Layer);
    RegisterPhysicalVolume(
        new G4PVPlacement(0, G4ThreeVector(-offsetLayer[ilay], 0, z_start + thicknessLayer[ilay] / 2), Log_Layer,
                          "physical_" + layer_name, log_sensor_stack, false, 0, overlapcheck_sector),
        layerActive[ilay]);
    z_start += thicknessLayer[ilay];
    m_DisplayAction->AddVolume(Log_Layer, "TTLLayers");
  }

  new G4PVReplica("TTL_Replicas_Modules", log_sensor_stack, log_sensor_ladder,
                  kYAxis, 6, baseplate_length);
  G4double offsety = diameter_coolingtube/2 + cooling_plate_height + thicknessDet / 2;
  G4double leftedge = sin(M_PI / 12.) * (rCenter + diameter_coolingtube/2 + cooling_plate_height );
  G4double leftedgebottom = sin(M_PI / 12.) * (rCenter - diameter_coolingtube/2 - cooling_plate_height);

  G4RotationMatrix *rotationSensor = new G4RotationMatrix();
  rotationSensor->rotateZ(-M_PI);
  // top side
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t1", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + baseplate_width / 2 + 2*baseSH_width_top + baseplate_width , 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t2", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + baseplate_width / 2 + 2*baseSH_width_top + 2* baseplate_width, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t3", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + baseplate_width / 2 + 4*baseSH_width_top + 3* baseplate_width, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t4", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + baseplate_width / 2 + 4*baseSH_width_top + 4* baseplate_width, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t5", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + baseplate_width / 2 + 6*baseSH_width_top + 5* baseplate_width, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t6", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + baseplate_width / 2 + 6*baseSH_width_top + 6* baseplate_width, 0, offsety), log_sensor_ladder,
                                            "physical_sensor_ladder_t7", log_module_envelope, false, 0, overlapcheck_sector),  false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + baseplate_width / 2 + 8*baseSH_width_top + 7* baseplate_width, 0, offsety), log_sensor_ladder,
                                            "physical_sensor_ladder_t8", log_module_envelope, false, 0, overlapcheck_sector),  false);

  // bottom side
  G4double offsetyDown = diameter_coolingtube/2 + cooling_plate_height + thicknessDet / 2;
  G4RotationMatrix *rotationSensorDown = new G4RotationMatrix();
  rotationSensorDown->rotateY(-M_PI);
  G4RotationMatrix *rotationSensorFlip = new G4RotationMatrix();
  rotationSensorFlip->rotateY(-M_PI);
  rotationSensorFlip->rotateZ(-M_PI);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedgebottom + baseSH_width + baseplate_width/2, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b1", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedgebottom + baseSH_width + baseplate_width/2 + baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b2", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedgebottom + 3*baseSH_width + baseplate_width/2 + 2*baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b3", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedgebottom + 3*baseSH_width + baseplate_width/2 + 3*baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b4", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedgebottom + 5*baseSH_width + baseplate_width/2 + 4*baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b5", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedgebottom + 5*baseSH_width + baseplate_width/2 + 5*baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b6", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedgebottom + 7*baseSH_width + baseplate_width/2 + 6*baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                            "physical_sensor_ladder_b7", log_module_envelope, false, 0, overlapcheck_sector),  false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedgebottom + 7*baseSH_width + baseplate_width/2 + 7*baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                            "physical_sensor_ladder_b8", log_module_envelope, false, 0, overlapcheck_sector),  false);
  // SERVICE HYBRID
  const int nLayers_SH = 4;
  std::string strLayerName_SH[nLayers_SH] = {
      "ThermalPad",
      "HighSpeedBoard",
      "ConnectorSpace",
      "Powerboard"};
  G4Material *materialLayer_SH[nLayers_SH] = {
      GetDetectorMaterial("G4_GRAPHITE"),
      GetDetectorMaterial("G4_POLYSTYRENE"),
      GetDetectorMaterial("G4_AIR"),
      GetDetectorMaterial("G4_POLYSTYRENE")};
  G4double thicknessLayer_SH[nLayers_SH] = {
      0.25 * mm,
      1.00 * mm,
      1.50 * mm,
      3.10 * mm};
  G4double widthLayer_SH[nLayers_SH] = {
      baseSH_width - 0.2 * mm,
      baseSH_width - 0.2 * mm,
      baseSH_width - 0.35 * mm,
      baseSH_width};
  G4double offsetLayer_SH[nLayers_SH] = {
      0.2 * mm / 2,
      0.2 * mm / 2,
      0.35 * mm / 2,
      0};
  G4double lengthLayer_SH[nLayers_SH] = {
      baseplate_length,
      baseplate_length,
      baseplate_length,
      baseplate_length};
  bool layerActive_SH[nLayers_SH] = {
      false, false, false, false};

  G4double thicknessDet_SH = 0;
  for (int ilay = 0; ilay < nLayers_SH; ilay++)
  {
    thicknessDet_SH += thicknessLayer_SH[ilay];
  }

  // Sensor Ladder (6 Sensors)
  G4VSolid *sol_SH_ladder = new G4Box("sol_SH_ladder",
                                      baseSH_width / 2,
                                      segmentlength / 2,
                                      thicknessDet_SH / 2);
  G4LogicalVolume *log_SH_ladder = new G4LogicalVolume(sol_SH_ladder, Air, "log_SH_ladder");
  m_DisplayAction->AddVolume(log_SH_ladder, "SHLadder");

  G4VSolid *sol_SH_stack = new G4Box("sol_SH_stack",
                                      baseSH_width / 2,
                                      baseplate_length / 2,
                                      thicknessDet_SH / 2);
  G4LogicalVolume *log_SH_stack = new G4LogicalVolume(sol_SH_stack, Air, "log_SH_stack");
  m_DisplayAction->AddVolume(log_SH_stack, "SHStack");

  double z_start_SH = -thicknessDet_SH / 2;
  for (int ilay = 0; ilay < nLayers_SH; ilay++)
  {
    const std::string layer_name = "SH_stack_" + strLayerName_SH[ilay];
    const std::string layer_name_Solid = "sol_" + layer_name;

    G4VSolid *sol_Module_Layer_Raw = new G4Box(layer_name_Solid + "_Raw",
                                                widthLayer_SH[ilay] / 2,
                                                lengthLayer_SH[ilay] / 2,
                                                thicknessLayer_SH[ilay] / 2);

    G4LogicalVolume *Log_Layer = new G4LogicalVolume(sol_Module_Layer_Raw,materialLayer_SH[ilay], layer_name + "_Log");
    RegisterLogicalVolume(Log_Layer);
    RegisterPhysicalVolume(
        new G4PVPlacement(0, G4ThreeVector(-offsetLayer_SH[ilay], 0, z_start_SH + thicknessLayer_SH[ilay] / 2), Log_Layer,
                          "physical_" + layer_name, log_SH_stack, false, 0, overlapcheck_sector),
        layerActive_SH[ilay]);
    z_start_SH += thicknessLayer_SH[ilay];
    m_DisplayAction->AddVolume(Log_Layer, "SHLayers");
  }

  new G4PVReplica("SH_Replicas_Modules", log_SH_stack, log_SH_ladder,
                  kYAxis, 6, baseplate_length);

  G4double offsety_SH = diameter_coolingtube/2 + cooling_plate_height + thicknessDet_SH / 2;
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 1 * baseplate_width + baseSH_width_top - baseSH_width/2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t1", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 1 * baseplate_width + 1 * baseSH_width_top + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t2", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 3 * baseplate_width + 3 * baseSH_width_top - baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t3", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 3 * baseplate_width + 3 * baseSH_width_top + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t4", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 5 * baseplate_width + 5 * baseSH_width_top - baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t5", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 5 * baseplate_width + 5 * baseSH_width_top + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t6", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 7 * baseplate_width + 7 * baseSH_width_top - baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                            "physical_SH_ladder_t7", log_module_envelope, false, 0, overlapcheck_sector),  false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 7 * baseplate_width + 7 * baseSH_width_top + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                            "physical_SH_ladder_t8", log_module_envelope, false, 0, overlapcheck_sector),  false);

  G4double offsetyDown_SH = diameter_coolingtube/2 + cooling_plate_height + thicknessDet_SH / 2;
  // RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedgebottom + 0 * baseplate_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
  //                                          "physical_SH_ladder_b1", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedgebottom + 2 * baseplate_width + 1 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b2", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedgebottom + 2 * baseplate_width + 2 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b3", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedgebottom + 4 * baseplate_width + 3 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b4", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedgebottom + 4 * baseplate_width + 4 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b5", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedgebottom + 6 * baseplate_width + 5 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b6", log_module_envelope, false, 0, overlapcheck_sector), false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedgebottom + 6 * baseplate_width + 6 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                            "physical_SH_ladder_b7", log_module_envelope, false, 0, overlapcheck_sector),  false);

  // RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedgebottom + 8 * baseplate_width + 7 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
  //                                           "physical_SH_ladder_b8", log_module_envelope, false, 0, overlapcheck_sector), false);


  RegisterLogicalVolume(log_module_envelope);
  m_DisplayAction->AddVolume(log_module_envelope, "ModuleEnvelope");
  G4double moduleShift = 3 * mm; //to avoid overlaps

  for (int isec = 0; isec < 12; isec++)
  {
    // if(isec!=3 && isec!=4)continue; // NOTE REMOVE
    // if(isec!=3)continue; // NOTE REMOVE
    G4RotationMatrix *motherrot = new G4RotationMatrix();
    motherrot->rotateX(M_PI / 2);
    motherrot->rotateY(M_PI);
    motherrot->rotateZ(M_PI + (isec - 3) * 2 * M_PI / 12.);
    // // central segments
    G4ThreeVector vec_central_transl((rCenter*cos(M_PI / 12.)+moduleShift) * cos(isec * 2 * M_PI / 12.), (rCenter*cos(M_PI / 12.)+moduleShift) * sin(isec * 2 * M_PI / 12.), place_z);
    assemblyDetector->AddPlacedVolume( log_module_envelope,vec_central_transl,motherrot);
    for (int ilen = 1; ilen < ((detlength / 2 - segmentlength / 2) / segmentlength); ilen++)
    {
      // forward segments
      G4ThreeVector vec_det_fwdlayers_transl((rCenter*cos(M_PI / 12.)+moduleShift) * cos(isec * 2 * M_PI / 12.), (rCenter*cos(M_PI / 12.)+moduleShift) * sin(isec * 2 * M_PI / 12.), place_z + ilen * segmentlength);
      assemblyDetector->AddPlacedVolume(log_module_envelope, vec_det_fwdlayers_transl, motherrot);
      // backward segments
      G4ThreeVector vec_det_bcklayers_transl((rCenter*cos(M_PI / 12.)+moduleShift) * cos(isec * 2 * M_PI / 12.), (rCenter*cos(M_PI / 12.)+moduleShift) * sin(isec * 2 * M_PI / 12.), place_z -ilen * segmentlength);
      assemblyDetector->AddPlacedVolume(log_module_envelope, vec_det_bcklayers_transl, motherrot);


    }


    // std::cout << "x : " << (rCenter*cos(M_PI / 12.)+moduleShift) * cos(isec * 2 * M_PI / 12.) << "\ty: " << (rCenter*cos(M_PI / 12.)+moduleShift) * sin(isec * 2 * M_PI / 12.) << "\tdistance: " << sqrt(pow((rCenter*cos(M_PI / 12.)+moduleShift) * cos(isec * 2 * M_PI / 12.),2)+pow((rCenter*cos(M_PI / 12.)+moduleShift) * sin(isec * 2 * M_PI / 12.),2)) << std::endl;

    // RegisterPhysicalVolume(new G4PVPlacement(motherrot, G4ThreeVector((rCenter*cos(M_PI / 12.)+moduleShift) * cos(isec * 2 * M_PI / 12.), (rCenter*cos(M_PI / 12.)+moduleShift) * sin(isec * 2 * M_PI / 12.), place_z), log_module_envelope,   "Mother_Segment_Raw_Physical_Center_" + std::to_string(isec), logicWorld, false, 0, overlapcheck_sector),false);
    // for (int ilen = 1; ilen < ((detlength / 2 - segmentlength / 2) / segmentlength); ilen++)
    // {
    //   // // forward segments
    //   RegisterPhysicalVolume(new G4PVPlacement(motherrot, G4ThreeVector((rCenter*cos(M_PI / 12.)+moduleShift) * cos(isec * 2 * M_PI / 12.), (rCenter*cos(M_PI / 12.)+moduleShift) * sin(isec * 2 * M_PI / 12.), place_z + ilen * segmentlength), log_module_envelope, "Mother_Segment_Raw_Physical_Fwd_" + std::to_string(isec) + "_" + std::to_string(ilen), logicWorld, false, 0, overlapcheck_sector),  false);
    //   // // backward segments
    //   RegisterPhysicalVolume(new G4PVPlacement(motherrot, G4ThreeVector((rCenter*cos(M_PI / 12.)+moduleShift) * cos(isec * 2 * M_PI / 12.), (rCenter*cos(M_PI / 12.)+moduleShift) * sin(isec * 2 * M_PI / 12.), place_z -ilen * segmentlength), log_module_envelope,"Mother_Segment_Raw_Physical_Bwd_" + std::to_string(isec) + "_" + std::to_string(ilen), logicWorld, false, 0, overlapcheck_sector), false);
    // }
  }
  assemblyDetector->MakeImprint( logicWorld, detzvec,0 );

  // RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, place_z), DetectorLog_Det,
  //                                          name_base + "_Physical", logicWorld, false, 0, overlapcheck_sector));
}

void PHG4TTLDetector::BuildForwardTTL(G4LogicalVolume *logicWorld)
{
  // MATERIALS
  G4double density;
  G4double a;
  G4int z;
  G4int ncomponents;
  G4int natoms;
  G4String symbol;

  m_SteppingAction->SetNPhiModules(0);
  m_SteppingAction->SetIsForwardTTL(true);


  G4Element *elH = new G4Element("Hydrogen", symbol = "H", 1., 1.01 * g / mole);
  G4Element *elC = new G4Element("Carbon", symbol = "C", 6., 12.01 * g / mole);
  G4Element *elN = new G4Element("Nitrogen", symbol = "N", 7., 14.01 * g / mole);
  G4Element *elO = new G4Element("Oxygen", symbol = "O", 8., 16.00 * g / mole);
  G4Material *mat_Epoxy = GetDetectorMaterial("EpoxyTTL",false);
  if (!mat_Epoxy)
  {
    mat_Epoxy = new G4Material("EpoxyTTL", density = 1.16 * g / cm3, natoms = 4);
    mat_Epoxy->AddElement(elH, 32);  // Hydrogen
    mat_Epoxy->AddElement(elN, 2);   // Nitrogen
    mat_Epoxy->AddElement(elO, 4);   // Oxygen
    mat_Epoxy->AddElement(elC, 15);  // Carbon
    // G4Material *mat_Epoxy = GetDetectorMaterial("Epoxy");
  }
  G4Material *mat_ALN = GetDetectorMaterial("AluminiumNitrate",false);
  if (!mat_ALN)
  {
    mat_ALN = new G4Material("AluminiumNitrate", density = 3.255 * g / cm3, ncomponents = 2);
    // G4Material *mat_ALN = new G4Material("AluminiumNitrate", density = 3.255 * g / cm3, ncomponents = 2);
    mat_ALN->AddElement(G4Element::GetElement("Al"), 1);
    mat_ALN->AddElement(G4Element::GetElement("N"), 1);
  }
  G4Material *mat_Solder_Tin = GetDetectorMaterial("Tin",false);
  if (!mat_Solder_Tin)
  {
    mat_Solder_Tin = new G4Material("Tin", z = 50., a = 118.7 * g / mole, density = 7.310 * g / cm3);
  }

  G4double det_height = 2.0 * cm;

  //Create the envelope = 'world volume' for the calorimeter
  G4Material *Air = GetDetectorMaterial("G4_AIR");

  G4double rMin = m_Params->get_double_param("rMin");
  G4double rMax = m_Params->get_double_param("rMax");
  G4double xoffset = m_Params->get_double_param("offset_x");
  // G4double xoffset = 6*cm;

  G4VSolid *beampipe_cutout = new G4Cons("ttl_beampipe_cutout",
                                         0, rMin,
                                         0, rMin,
                                         det_height * 2 / 2.0,
                                         0, 2 * M_PI);
  G4VSolid *ttl_envelope_solid = new G4Cons("ttl_envelope_solid_cutout",
                                            0, rMax,
                                            0, rMax,
                                            det_height / 2.0,
                                            0, 2 * M_PI);
  ttl_envelope_solid = new G4SubtractionSolid(G4String("ttl_envelope_solid"), ttl_envelope_solid, beampipe_cutout, 0, G4ThreeVector(xoffset, 0, 0.));

  const G4Transform3D transform_Det_to_Hall =
      G4RotateX3D(-m_Params->get_double_param("polar_angle")) * G4TranslateZ3D(
                                                                    m_Params->get_double_param("place_z") + det_height / 2);

  m_SteppingAction->SetZPositionFwd(m_Params->get_double_param("place_z") + det_height / 2);
  G4LogicalVolume *DetectorLog_Det = new G4LogicalVolume(ttl_envelope_solid, Air, name_base + "_Log");
  RegisterLogicalVolume(DetectorLog_Det);

  RegisterPhysicalVolume(new G4PVPlacement(G4RotateZ3D(2 * pi) * transform_Det_to_Hall, DetectorLog_Det,
                                           name_base + "_Physical", logicWorld, false, 0, overlapcheck_sector));
  m_DisplayAction->AddVolume(DetectorLog_Det, "DetectorBoxFwd");

  G4double cooling_plate_height = 1 * mm;
  if(m_Params->get_double_param("cooling_plate_height")>0) cooling_plate_height = m_Params->get_double_param("cooling_plate_height");

  G4double diameter_coolingtube = 5*mm;
  G4double wallthickness_coolingtube = 1 * mm;

  G4VSolid *sol_module_envelope = new G4Cons("sol_module_envelope_cutout",
                                          0, rMax,
                                          0, rMax,
                                          det_height / 2.0,
                                          0, 2 * M_PI);
  sol_module_envelope = new G4SubtractionSolid(G4String("sol_module_envelope_"), sol_module_envelope, beampipe_cutout, 0, G4ThreeVector(xoffset, 0, 0.));
  G4LogicalVolume *log_module_envelope = new G4LogicalVolume(sol_module_envelope, Air, "log_module_envelope");
  m_DisplayAction->AddVolume(log_module_envelope, "ModuleEnvelope");
  RegisterLogicalVolume(log_module_envelope);
  // G4RotationMatrix *rot1 = new G4RotationMatrix();
  // rot1->rotateZ(0.5 * isec * M_PI);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, 0), log_module_envelope,
                                            "physical_module_envelope2", DetectorLog_Det, false, 0, overlapcheck_sector));

  G4VSolid *sol_cooling_plate = new G4Cons("sol_cooling_plate_cutout",
                                            0, rMax,
                                            0, rMax,
                                            cooling_plate_height / 2.0,
                                            0, 2 * M_PI);
  sol_cooling_plate = new G4SubtractionSolid(G4String("sol_cooling_plate"), sol_cooling_plate, beampipe_cutout, 0, G4ThreeVector(xoffset, 0, 0.));

  G4LogicalVolume *log_cooling_plate = new G4LogicalVolume(sol_cooling_plate, GetDetectorMaterial("G4_Al"), "log_cooling_plate");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, diameter_coolingtube/2+cooling_plate_height/2), log_cooling_plate,
                                            "physical_cooling_plate_f", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, -diameter_coolingtube/2-cooling_plate_height/2), log_cooling_plate,
                                            "physical_cooling_plate_b", log_module_envelope, false, 0, overlapcheck_sector), false);
  RegisterLogicalVolume(log_cooling_plate);
  m_DisplayAction->AddVolume(log_cooling_plate, "CoolingPlate");


  // Sensor Module:
  G4double sensor_width = 21.2 * mm;
  G4double sensor_length = 42.0 * mm;
  G4double baseplate_length = 43.1 * mm;
  G4double baseplate_width = 56.5 * mm / 2;

  const int nLayers = 8;
  std::string strLayerName[nLayers] = {
      "ThermalPad",
      "ALN",
      "LairdFilm",
      "ROC",
      "Solder",
      "Sensor",
      "Epoxy",
      "AIN"};
  G4Material *materialLayer[nLayers] = {GetDetectorMaterial("G4_GRAPHITE"), mat_ALN, GetDetectorMaterial("G4_GRAPHITE"), GetDetectorMaterial("G4_PLEXIGLASS"), mat_Solder_Tin, GetDetectorMaterial("G4_Si"), mat_Epoxy, mat_ALN};
  G4double thicknessLayer[nLayers] = {
      0.25 * mm, 0.79 * mm, 0.08 * mm, 0.25 * mm, 0.03 * mm, 0.3 * mm, 0.08 * mm, 0.51 * mm};
  G4double widthLayer[nLayers] = {
      baseplate_width - 0.3 * mm,
      baseplate_width,
      sensor_width + 1 * mm,
      sensor_width + 1 * mm,
      sensor_width - 0.2 * mm,
      sensor_width,
      sensor_width,
      baseplate_width - 4 * mm};
  G4double offsetLayer[nLayers] = {
      0.3 * mm / 2,
      0,
      (baseplate_width - widthLayer[2]) / 2 - 0.2 * mm,
      (baseplate_width - widthLayer[3]) / 2 - 0.2 * mm,
      (baseplate_width - widthLayer[4]) / 2 - 0.1 * mm,
      (baseplate_width - widthLayer[5] - 0.1 * mm) / 2,
      (baseplate_width - widthLayer[6] - 0.1 * mm) / 2,
      4 * mm / 2 - 0.2 * mm};
  G4double lengthLayer[nLayers] = {
      baseplate_length - 0.2 * mm,
      baseplate_length,
      sensor_length + 0.2 * mm,
      sensor_length + 0.2 * mm,
      sensor_length - 0.2 * mm,
      sensor_length,
      sensor_length,
      baseplate_length - 0.2 * mm};
  bool layerActive[nLayers] = {
      false, false, false, false, false, true, false, false};

  G4double thicknessDet = 0;
  for (int ilay = 0; ilay < nLayers; ilay++)
  {
    thicknessDet += thicknessLayer[ilay];
  }

  G4double segmentlength = baseplate_length;  //(detlength - 10 * cm) / 6;//
  G4VSolid *sol_sensor_stack = new G4Box("sol_sensor_stack",
                                         baseplate_length / 2,
                                         baseplate_width / 2,
                                         thicknessDet / 2);
  G4LogicalVolume *log_sensor_stack = new G4LogicalVolume(sol_sensor_stack, Air, "log_sensor_stack");
  m_DisplayAction->AddVolume(log_sensor_stack, "SensorStack");

  double z_start = -thicknessDet / 2;
  for (int ilay = 0; ilay < nLayers; ilay++)
  {
    const std::string layer_name = "sensor_stack_" + strLayerName[ilay];
    const std::string layer_name_Solid = "sol_" + layer_name;

    G4VSolid *sol_Module_Layer_Raw = new G4Box(layer_name_Solid + "_Raw",
                                               lengthLayer[ilay] / 2,
                                               widthLayer[ilay] / 2,
                                               thicknessLayer[ilay] / 2);

    G4LogicalVolume *Log_Layer = new G4LogicalVolume(sol_Module_Layer_Raw,  //
                                                     materialLayer[ilay], layer_name + "_Log");
    RegisterLogicalVolume(Log_Layer);
    RegisterPhysicalVolume(
        new G4PVPlacement(0, G4ThreeVector(0, -offsetLayer[ilay], z_start + thicknessLayer[ilay] / 2), Log_Layer,
                          "physical_" + layer_name, log_sensor_stack, false, 0, overlapcheck_sector),
        layerActive[ilay]);
    z_start += thicknessLayer[ilay];
    m_DisplayAction->AddVolume(Log_Layer, "TTLLayers");
  }

  // SERVICE HYBRID
  G4double baseSH_width = baseplate_width / 2;
  const int nLayers_SH = 4;
  std::string strLayerName_SH[nLayers_SH] = {
      "ThermalPad",
      "HighSpeedBoard",
      "ConnectorSpace",
      "Powerboard"};
  G4Material *materialLayer_SH[nLayers_SH] = {
      GetDetectorMaterial("G4_GRAPHITE"),
      GetDetectorMaterial("G4_POLYSTYRENE"),
      GetDetectorMaterial("G4_AIR"),
      GetDetectorMaterial("G4_POLYSTYRENE")};
  G4double thicknessLayer_SH[nLayers_SH] = {
      0.25 * mm,
      1.00 * mm,
      1.50 * mm,
      3.10 * mm};
  G4double widthLayer_SH[nLayers_SH] = {
      baseSH_width - 0.2 * mm,
      baseSH_width - 0.2 * mm,
      baseSH_width - 0.35 * mm,
      baseSH_width};
  G4double offsetLayer_SH[nLayers_SH] = {
      0.2 * mm / 2,
      0.2 * mm / 2,
      0.35 * mm / 2,
      0};
  G4double lengthLayer_SH[nLayers_SH] = {
      baseplate_length,
      baseplate_length,
      baseplate_length,
      baseplate_length};
  bool layerActive_SH[nLayers_SH] = {
      false, false, false, false};

  G4double thicknessDet_SH = 0;
  for (int ilay = 0; ilay < nLayers_SH; ilay++)
  {
    thicknessDet_SH += thicknessLayer_SH[ilay];
  }

  G4VSolid *sol_SH_stack = new G4Box("sol_SH_stack",
                                     baseplate_length / 2,
                                     baseSH_width / 2,
                                     thicknessDet_SH / 2);
  G4LogicalVolume *log_SH_stack = new G4LogicalVolume(sol_SH_stack, Air, "log_SH_stack");
  m_DisplayAction->AddVolume(log_SH_stack, "SHStack");

  double z_start_SH = -thicknessDet_SH / 2;
  for (int ilay = 0; ilay < nLayers_SH; ilay++)
  {
    const std::string layer_name = "SH_stack_" + strLayerName_SH[ilay];
    const std::string layer_name_Solid = "sol_" + layer_name;

    G4VSolid *sol_Module_Layer_Raw = new G4Box(layer_name_Solid + "_Raw",
                                               lengthLayer_SH[ilay] / 2,
                                               widthLayer_SH[ilay] / 2,
                                               thicknessLayer_SH[ilay] / 2);

    G4LogicalVolume *Log_Layer = new G4LogicalVolume(sol_Module_Layer_Raw,  //
                                                     materialLayer_SH[ilay], layer_name + "_Log");
    RegisterLogicalVolume(Log_Layer);
    RegisterPhysicalVolume(
        new G4PVPlacement(0, G4ThreeVector(0,-offsetLayer_SH[ilay], z_start_SH + thicknessLayer_SH[ilay] / 2), Log_Layer,
                          "physical_" + layer_name, log_SH_stack, false, 0, overlapcheck_sector),
        layerActive_SH[ilay]);
    z_start_SH += thicknessLayer_SH[ilay];
    m_DisplayAction->AddVolume(Log_Layer, "SHLayers");
  }

  G4double fullsensor_width = baseSH_width+baseplate_width;
  G4VSolid *sol_sensor_and_readout = new G4Box("sol_sensor_and_readout",
                                          segmentlength / 2,
                                          fullsensor_width / 2,
                                          thicknessDet_SH / 2);
  G4LogicalVolume *log_sensor_and_readout = new G4LogicalVolume(sol_sensor_and_readout, Air, "log_sensor_and_readout");
  m_DisplayAction->AddVolume(log_sensor_and_readout, "SensorAndReadoutLadder");

  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet/2),
                      log_sensor_stack, "SensorPlacedPhysical", log_sensor_and_readout, false, 0, overlapcheck_sector),false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, fullsensor_width/2 - baseplate_width/4, 0),
                      log_SH_stack, "ServiceHybridPlacedPhysical", log_sensor_and_readout, false, 0, overlapcheck_sector),false);

  G4double offsetzFront = diameter_coolingtube/2 + cooling_plate_height + thicknessDet_SH / 2;
  G4double offsetzBack = -diameter_coolingtube/2 - cooling_plate_height - thicknessDet_SH / 2;
  // number of towers in radial direction (on y axis)
  int rowYdir = (int) ( (rMax-(fullsensor_width/2)) / fullsensor_width);
  for(int row=rowYdir;row>=-rowYdir;row--){
  // for(int row=0;row>=-rowYdir;row--){
    // pythagoras -> get available length in circular mother volume for towers
    // divide given length by tower width -> get number of towers that can be placed
    int numSensorsRow = (int) ( ( 2* sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width + fullsensor_width/2) ,2)) ) / segmentlength );
    if(numSensorsRow==0) continue;
    // we want an odd number of towers to be symmetrically centered around 0
    if ( numSensorsRow % 2 == 0) numSensorsRow-=1;

    if( ( (abs(row)*fullsensor_width) -(fullsensor_width/2)) < rMin ){
      if(xoffset!=0){
        // pythagoras -> get available length in circular mother volume for towers
        // divide given length by tower width -> get number of towers that can be placed
        // int numSensorsRowInner = ceil( ( 2* sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) ) / segmentlength );
        int numSensorLeftAdd = ceil( (xoffset -(segmentlength/2) - sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) ) / segmentlength );
        int numSensorRightAdd = ceil( (xoffset -(segmentlength/2) + sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) ) / segmentlength );
        // we want an odd number of towers to be symmetrically centered around 0
        // create mother volume with space for numSensorsRow towers along x-axis
        auto TTLDetRowLeftSolid    = new G4Box("TTLDetRowLeftBox" + std::to_string(row), (((numSensorsRow-1) /2 + numSensorLeftAdd)) * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        auto TTLDetRowLeftLogical  = new G4LogicalVolume(TTLDetRowLeftSolid,Air,"TTLDetRowLeftLogical" + std::to_string(row));
        // replicate singletower tower design numSensorsRow times along x-axis
        new G4PVReplica("TTLDetRowLeftPhysical" + std::to_string(row),log_sensor_and_readout,TTLDetRowLeftLogical,
                        kXAxis,((numSensorsRow-1) /2 + numSensorLeftAdd),segmentlength);
        m_DisplayAction->AddVolume(TTLDetRowLeftLogical, "StripBox");

        G4RotationMatrix *rotationSensor = new G4RotationMatrix();
        rotationSensor->rotateZ(-M_PI);
        new G4PVPlacement(row%2==0 ? rotationSensor : 0, G4ThreeVector( - ( (((numSensorsRow-1) /2 + numSensorLeftAdd)) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), (row*fullsensor_width), offsetzFront),
                      TTLDetRowLeftLogical, "TTLDetRowLeftPlacedFront" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector);

        G4RotationMatrix *rotationSensorBck1 = new G4RotationMatrix();
        rotationSensorBck1->rotateX(-M_PI);
        G4RotationMatrix *rotationSensorBck2 = new G4RotationMatrix();
        rotationSensorBck2->rotateZ(-M_PI);
        rotationSensorBck2->rotateX(-M_PI);
        new G4PVPlacement(row%2==0 ? rotationSensorBck2 : rotationSensorBck1, G4ThreeVector( - ( (((numSensorsRow-1) /2 + numSensorLeftAdd)) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), (row*fullsensor_width), offsetzBack),
                      TTLDetRowLeftLogical, "TTLDetRowLeftPlacedBack" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector);

        G4VSolid *sol_cutout_tube_left = new G4Box("sol_cutout_tube_left" + std::to_string(row),
                                                1.3*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) + xoffset)/2,
                                                (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                                (diameter_coolingtube - 2*wallthickness_coolingtube) / 2);
        G4VSolid *sol_cooling_tube_left = new G4Box("sol_cooling_tube_left_tmp" + std::to_string(row),
                                                0.86*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) + xoffset)/2,
                                                diameter_coolingtube / 2,
                                                diameter_coolingtube / 2);
        sol_cooling_tube_left = new G4SubtractionSolid(G4String("sol_cooling_tube_left" + std::to_string(row)), sol_cooling_tube_left, sol_cutout_tube_left, 0, G4ThreeVector(0,0,0));
        G4LogicalVolume *Log_cooling_tube_left = new G4LogicalVolume(sol_cooling_tube_left,  //
                                                              GetDetectorMaterial("G4_Al"), "Log_cooling_tube_left" + std::to_string(row));
        RegisterLogicalVolume(Log_cooling_tube_left);
        m_DisplayAction->AddVolume(Log_cooling_tube_left, "Cooling_tube");

        G4VSolid *sol_water_cooling_left = new G4Box("sol_water_cooling_left" + std::to_string(row),
                                                0.85*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) + xoffset)/2,
                                                0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                                0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2);
        G4LogicalVolume *Log_water_cooling_left = new G4LogicalVolume(sol_water_cooling_left,  //
                                                              GetDetectorMaterial("G4_WATER"), "Log_water_cooling_left" + std::to_string(row));
        RegisterLogicalVolume(Log_water_cooling_left);
        m_DisplayAction->AddVolume(Log_water_cooling_left, "Water_cooling");

        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(- ( (((numSensorsRow-1) /2 + numSensorLeftAdd)) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_cooling_tube_left,
                                        "cooling_tube_left_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);
        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(- ( (((numSensorsRow-1) /2 + numSensorLeftAdd)) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_water_cooling_left,
                                        "cooling_water_left_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);


        // // create mother volume with space for numSensorsRow towers along x-axis
        auto TTLDetRowRightSolid    = new G4Box("TTLDetRowRightBox" + std::to_string(row), (((numSensorsRow-1) /2 - numSensorRightAdd)) * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        auto TTLDetRowRightLogical  = new G4LogicalVolume(TTLDetRowRightSolid,Air,"TTLDetRowRightLogical" + std::to_string(row));
        m_DisplayAction->AddVolume(TTLDetRowRightLogical, "StripBox");
        // // replicate singletower tower design numSensorsRow times along x-axis
        new G4PVReplica("TTLDetRowRightPhysical" + std::to_string(row),log_sensor_and_readout,TTLDetRowRightLogical,
                        kXAxis,((numSensorsRow-1) /2 - numSensorRightAdd ),segmentlength);

        new G4PVPlacement(row%2==0 ? rotationSensor : 0, G4ThreeVector(( (((numSensorsRow-1) /2 - numSensorRightAdd)) * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), (row*fullsensor_width), offsetzFront),
                      TTLDetRowRightLogical, "TTLDetRowRightPlacedFront" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector);

        new G4PVPlacement(row%2==0 ? rotationSensorBck2 : rotationSensorBck1, G4ThreeVector(( (((numSensorsRow-1) /2 - numSensorRightAdd)) * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), (row*fullsensor_width), offsetzBack),
                      TTLDetRowRightLogical, "TTLDetRowRightPlacedBack" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector);


        G4VSolid *sol_cutout_tube_right = new G4Box("sol_cutout_tube_right" + std::to_string(row),
                                                1.3*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) - xoffset)/2,
                                                (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                                (diameter_coolingtube - 2*wallthickness_coolingtube) / 2);
        G4VSolid *sol_cooling_tube_right = new G4Box("sol_cooling_tube_right_tmp" + std::to_string(row),
                                                0.86*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) - xoffset)/2,
                                                diameter_coolingtube / 2,
                                                diameter_coolingtube / 2);
        sol_cooling_tube_right = new G4SubtractionSolid(G4String("sol_cooling_tube_right" + std::to_string(row)), sol_cooling_tube_right, sol_cutout_tube_right, 0, G4ThreeVector(0,0,0));
        G4LogicalVolume *Log_cooling_tube_right = new G4LogicalVolume(sol_cooling_tube_right,  //
                                                              GetDetectorMaterial("G4_Al"), "Log_cooling_tube_right" + std::to_string(row));
        RegisterLogicalVolume(Log_cooling_tube_right);
        m_DisplayAction->AddVolume(Log_cooling_tube_right, "Cooling_tube");

        G4VSolid *sol_water_cooling_right = new G4Box("sol_water_cooling_right" + std::to_string(row),
                                                0.85*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) - xoffset)/2,
                                                0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                                0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2);
        G4LogicalVolume *Log_water_cooling_right = new G4LogicalVolume(sol_water_cooling_right,  //
                                                              GetDetectorMaterial("G4_WATER"), "Log_water_cooling_right" + std::to_string(row));
        RegisterLogicalVolume(Log_water_cooling_right);
        m_DisplayAction->AddVolume(Log_water_cooling_right, "Water_cooling");

        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(( (((numSensorsRow-1) /2 - numSensorRightAdd)) * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_cooling_tube_right,
                                        "cooling_tube_right_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);
        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(( (((numSensorsRow-1) /2 - numSensorRightAdd)) * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_water_cooling_right,
                                        "cooling_water_right_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);

      } else {
        // pythagoras -> get available length in circular mother volume for towers
        // divide given length by tower width -> get number of towers that can be placed
        int numSensorsInner = ceil( ( 2* sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2)) ) / segmentlength );
        // we want an odd number of towers to be symmetrically centered around 0
        if ( numSensorsInner % 2 == 0) numSensorsInner+=1;
        // numSensorsInner+=2;
        // create mother volume with space for numSensors towers along x-axis
        auto TTLDetRowSolid    = new G4Box("TTLDetRowBox" + std::to_string(row), (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        auto TTLDetRowLogical  = new G4LogicalVolume(TTLDetRowSolid,Air,"TTLDetRowLogical" + std::to_string(row));
        m_DisplayAction->AddVolume(TTLDetRowLogical, "StripBox");
        // replicate singletower tower design numSensors times along x-axis
        new G4PVReplica("TTLDetRowPhysical" + std::to_string(row),log_sensor_and_readout,TTLDetRowLogical,
                        kXAxis,(numSensorsRow - numSensorsInner) / 2,segmentlength);

        G4RotationMatrix *rotationSensor = new G4RotationMatrix();
        rotationSensor->rotateZ(-M_PI);
        RegisterPhysicalVolume(new G4PVPlacement(row%2==0 ? rotationSensor : 0, G4ThreeVector( - ( ( numSensorsInner / 2.0 ) * segmentlength ) - ( (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0 ), (row*fullsensor_width), offsetzFront),
                      TTLDetRowLogical, "TTLDetRowPhysicalPlacedFrontLeft_" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector),false);

        RegisterPhysicalVolume(new G4PVPlacement(row%2==0 ? rotationSensor : 0, G4ThreeVector( ( ( numSensorsInner / 2.0 ) * segmentlength ) + ( (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0 ), (row*fullsensor_width), offsetzFront),
                      TTLDetRowLogical, "TTLDetRowPhysicalPlacedFrontRight_" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector),false);

        G4RotationMatrix *rotationSensorBck1 = new G4RotationMatrix();
        rotationSensorBck1->rotateX(-M_PI);
        G4RotationMatrix *rotationSensorBck2 = new G4RotationMatrix();
        rotationSensorBck2->rotateZ(-M_PI);
        rotationSensorBck2->rotateX(-M_PI);

        RegisterPhysicalVolume(new G4PVPlacement(row%2==0 ? rotationSensorBck2 : rotationSensorBck1, G4ThreeVector( - ( ( numSensorsInner / 2.0 ) * segmentlength ) - ( (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0 ), (row*fullsensor_width), offsetzBack),
                      TTLDetRowLogical, "TTLDetRowPhysicalPlacedBckLeft_" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector),false);

        RegisterPhysicalVolume(new G4PVPlacement(row%2==0 ? rotationSensorBck2 : rotationSensorBck1, G4ThreeVector( ( ( numSensorsInner / 2.0 ) * segmentlength ) + ( (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0 ), (row*fullsensor_width), offsetzBack),
                      TTLDetRowLogical, "TTLDetRowPhysicalPlacedBckRight_" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector),false);

        // cooling
        G4VSolid *sol_cutout_tube_bothsides = new G4Box("sol_cutout_tube_bothsides" + std::to_string(row),
                                                1.3*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) )/2,
                                                (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                                (diameter_coolingtube - 2*wallthickness_coolingtube) / 2);
        G4VSolid *sol_cooling_tube_bothsides = new G4Box("sol_cooling_tube_bothsides_tmp" + std::to_string(row),
                                                0.96*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) )/2,
                                                diameter_coolingtube / 2,
                                                diameter_coolingtube / 2);
        sol_cooling_tube_bothsides = new G4SubtractionSolid(G4String("sol_cooling_tube_bothsides" + std::to_string(row)), sol_cooling_tube_bothsides, sol_cutout_tube_bothsides, 0, G4ThreeVector(0,0,0));
        G4LogicalVolume *Log_cooling_tube_bothsides = new G4LogicalVolume(sol_cooling_tube_bothsides,  //
                                                              GetDetectorMaterial("G4_Al"), "Log_cooling_tube_bothsides" + std::to_string(row));
        RegisterLogicalVolume(Log_cooling_tube_bothsides);
        m_DisplayAction->AddVolume(Log_cooling_tube_bothsides, "Cooling_tube");

        G4VSolid *sol_water_cooling_bothsides = new G4Box("sol_water_cooling_bothsides" + std::to_string(row),
                                                0.95*(sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)))/2,
                                                0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                                0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2);
        G4LogicalVolume *Log_water_cooling_bothsides = new G4LogicalVolume(sol_water_cooling_bothsides,  //
                                                              GetDetectorMaterial("G4_WATER"), "Log_water_cooling_bothsides" + std::to_string(row));
        RegisterLogicalVolume(Log_water_cooling_bothsides);
        m_DisplayAction->AddVolume(Log_water_cooling_bothsides, "Water_cooling");

        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(- ( ( numSensorsInner / 2.0 ) * segmentlength ) - ( (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0 ), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_cooling_tube_bothsides,
                                        "cooling_tube_bothsides_l_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);
        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(( ( numSensorsInner / 2.0 ) * segmentlength ) + ( (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0 ), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_cooling_tube_bothsides,
                                        "cooling_tube_bothsides_r_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);
        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(- ( ( numSensorsInner / 2.0 ) * segmentlength ) - ( (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0 ), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_water_cooling_bothsides,
                                        "cooling_water_bothsides_l_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);
        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(( ( numSensorsInner / 2.0 ) * segmentlength ) + ( (numSensorsRow - numSensorsInner) / 2 * segmentlength / 2.0 ), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_water_cooling_bothsides,
                                        "cooling_water_bothsides_r_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);
      }
    } else {
        // create mother volume with space for numSensorsRow towers along x-axis
        auto TTLDetRowSolid    = new G4Box("TTLDetRowBox" + std::to_string(row), numSensorsRow * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        auto TTLDetRowLogical  = new G4LogicalVolume(TTLDetRowSolid,Air,"TTLDetRowLogical" + std::to_string(row));
        m_DisplayAction->AddVolume(TTLDetRowLogical, "StripBox");
        // replicate singletower tower design numSensorsRow times along x-axis
        new G4PVReplica("TTLDetRowPhysicalReplica" + std::to_string(row),log_sensor_and_readout,TTLDetRowLogical,
                        kXAxis,numSensorsRow,segmentlength);

        G4RotationMatrix *rotationSensor = new G4RotationMatrix();
        rotationSensor->rotateZ(-M_PI);
        RegisterPhysicalVolume(new G4PVPlacement(row%2==0 ? rotationSensor : 0, G4ThreeVector(0, (row*fullsensor_width), offsetzFront),
                      TTLDetRowLogical, "TTLDetRowPhysicalFront" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),false);

        G4RotationMatrix *rotationSensorBck1 = new G4RotationMatrix();
        rotationSensorBck1->rotateX(-M_PI);
        G4RotationMatrix *rotationSensorBck2 = new G4RotationMatrix();
        rotationSensorBck2->rotateZ(-M_PI);
        rotationSensorBck2->rotateX(-M_PI);
        RegisterPhysicalVolume(new G4PVPlacement(row%2==0 ? rotationSensorBck2 : rotationSensorBck1, G4ThreeVector(0, (row*fullsensor_width), offsetzBack),
                      TTLDetRowLogical, "TTLDetRowPhysicalBack" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),false);


        G4VSolid *sol_cutout_tube = new G4Box("sol_cutout_tube" + std::to_string(row),
                                                2.5*sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) + diameter_coolingtube ,2))/2,
                                                (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                                (diameter_coolingtube - 2*wallthickness_coolingtube) / 2);
        G4VSolid *sol_cooling_tube = new G4Box("sol_cooling_tube_tmp" + std::to_string(row),
                                                0.98*2*sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) + diameter_coolingtube ,2))/2,
                                                diameter_coolingtube / 2,
                                                diameter_coolingtube / 2);
        sol_cooling_tube = new G4SubtractionSolid(G4String("sol_cooling_tube" + std::to_string(row)), sol_cooling_tube, sol_cutout_tube, 0, G4ThreeVector(0,0,0));
        G4LogicalVolume *Log_cooling_tube = new G4LogicalVolume(sol_cooling_tube,  //
                                                              GetDetectorMaterial("G4_Al"), "Log_cooling_tube" + std::to_string(row));
        RegisterLogicalVolume(Log_cooling_tube);
        m_DisplayAction->AddVolume(Log_cooling_tube, "Cooling_tube");

        G4VSolid *sol_water_cooling = new G4Box("sol_water_cooling" + std::to_string(row),
                                                0.97*2*sqrt(pow(rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) + diameter_coolingtube ,2))/2,
                                                0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
                                                0.99*(diameter_coolingtube - 2*wallthickness_coolingtube) / 2);
        G4LogicalVolume *Log_water_cooling = new G4LogicalVolume(sol_water_cooling,  //
                                                              GetDetectorMaterial("G4_WATER"), "Log_water_cooling" + std::to_string(row));
        RegisterLogicalVolume(Log_water_cooling);
        m_DisplayAction->AddVolume(Log_water_cooling, "Water_cooling");

        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(0, row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_cooling_tube,
                                        "cooling_tube_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);
        RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(0, row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), 0), Log_water_cooling,
                                        "cooling_water_Physical_" + std::to_string(row), log_module_envelope, false, 0, overlapcheck_sector),   false);
    }
  }
}


G4Material* PHG4TTLDetector::MakeCarbonFoamMaterial(){
  G4Material* carbon_foam = GetDetectorMaterial("C_FOAM_TTL", false);  // false suppresses warning that material does not exist
  if(!carbon_foam){
    G4double density;
    G4int ncomponents;
    carbon_foam = new G4Material("C_FOAM_TTL", density = 0.20 * g / cm3, ncomponents = 2); // VERY CONSERVATIVE DENSITY
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.97);
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 0.03);
  }
  return carbon_foam;

}

G4Material* PHG4TTLDetector::GetCarbonFiber()
{
  G4Material* carbonfiber = G4Material::GetMaterial("TTLCarbonFiber", false);  // false suppresses warning that material does not exist
  if (!carbonfiber)
  {
    G4double density_carbon_fiber = 1.44 * g / cm3;
    carbonfiber = new G4Material("TTLCarbonFiber", density_carbon_fiber, 1);
    carbonfiber->AddElement(G4Element::GetElement("C"), 1);
  }
  return carbonfiber;
}

G4LogicalVolume *
PHG4TTLDetector::RegisterLogicalVolume(G4LogicalVolume *v)
{
  if (!v)
  {
    std::cout
        << "PHG4TTLDetector::RegisterVolume - Error - invalid volume!"
        << std::endl;
    return v;
  }
  if (map_log_vol.find(v->GetName()) != map_log_vol.end())
  {
    std::cout << "PHG4TTLDetector::RegisterVolume - Warning - replacing " << v->GetName() << std::endl;
  }

  map_log_vol[v->GetName()] = v;

  return v;
}

G4PVPlacement *
PHG4TTLDetector::RegisterPhysicalVolume(G4PVPlacement *v,
                                        const bool active)
{
  if (!v)
  {
    std::cout
        << "PHG4TTLDetector::RegisterPhysicalVolume - Error - invalid volume!"
        << std::endl;
    return v;
  }

  phy_vol_idx_t id(v->GetName(), v->GetCopyNo());

  if (map_phy_vol.find(id) != map_phy_vol.end())
  {
    std::cout
        << "PHG4TTLDetector::RegisterPhysicalVolume - Warning - replacing "
        << v->GetName() << "[" << v->GetCopyNo() << "]" << std::endl;
  }

  map_phy_vol[id] = v;

  if (active)
    map_active_phy_vol[id] = v;

  return v;
}
