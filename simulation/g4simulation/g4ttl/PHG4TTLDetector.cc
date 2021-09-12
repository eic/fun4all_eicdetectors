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
    BuildBarrelTTL(logicWorld);
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

  G4Element *elH = new G4Element("Hydrogen", symbol = "H", 1., 1.01 * g / mole);
  G4Element *elC = new G4Element("Carbon", symbol = "C", 6., 12.01 * g / mole);
  G4Element *elN = new G4Element("Nitrogen", symbol = "N", 7., 14.01 * g / mole);
  G4Element *elO = new G4Element("Oxygen", symbol = "O", 8., 16.00 * g / mole);
  G4Material *mat_Epoxy = G4Material::GetMaterial("EpoxyTTL");
  if (!mat_Epoxy)
  {
    mat_Epoxy = new G4Material("EpoxyTTL", density = 1.16 * g / cm3, natoms = 4);
    mat_Epoxy->AddElement(elH, 32);  // Hydrogen
    mat_Epoxy->AddElement(elN, 2);   // Nitrogen
    mat_Epoxy->AddElement(elO, 4);   // Oxygen
    mat_Epoxy->AddElement(elC, 15);  // Carbon
    // G4Material *mat_Epoxy = G4Material::GetMaterial("Epoxy");
  }
  G4Material *mat_ALN = G4Material::GetMaterial("AluminiumNitrate");
  if (!mat_ALN)
  {
    mat_ALN = new G4Material("AluminiumNitrate", density = 3.255 * g / cm3, ncomponents = 2);
    // G4Material *mat_ALN = new G4Material("AluminiumNitrate", density = 3.255 * g / cm3, ncomponents = 2);
    mat_ALN->AddElement(G4Element::GetElement("Al"), 1);
    mat_ALN->AddElement(G4Element::GetElement("N"), 1);
  }
  G4Material *mat_Solder_Tin = G4Material::GetMaterial("Tin");
  if (!mat_Solder_Tin)
  {
    mat_Solder_Tin = new G4Material("Tin", z = 50., a = 118.7 * g / mole, density = 7.310 * g / cm3);
  }
  G4Material *Air = G4Material::GetMaterial("G4_AIR");

  // positions
  G4double rMin = m_Params->get_double_param("rMin");  // center location of Al support plate
  G4double det_height = 2.0 * cm;
  G4double place_z = m_Params->get_double_param("place_z");
  G4double detlength = m_Params->get_double_param("length");

  //Create the envelope = 'world volume' for the calorimeter
  G4VSolid *ttl_envelope_solid = new G4Cons("ttl_envelope_solid",
                                            rMin - det_height / 2 - 2 * cm, rMin + det_height / 2 + 2 * cm,
                                            rMin - det_height / 2 - 2 * cm, rMin + det_height / 2 + 2 * cm,
                                            detlength / 2.0,
                                            0, 2 * M_PI);

  G4LogicalVolume *DetectorLog_Det = new G4LogicalVolume(ttl_envelope_solid, Air, name_base + "_Log");
  RegisterLogicalVolume(DetectorLog_Det);
  m_DisplayAction->AddVolume(DetectorLog_Det, "FullEnvelope");

  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, place_z), DetectorLog_Det,
                                           name_base + "_Physical", logicWorld, false, 0, overlapcheck_sector));

  // Single module with length based on readout (contains 14 LGADs [counting across both sides] in x-direction and 6 in z-direction)
  G4double baseplate_length = 43.1 * mm;
  G4double baseplate_width = 56.5 * mm / 2;
  G4double segmentlength = 6 * baseplate_length;  //(detlength - 10 * cm) / 6;//m_Params->get_double_param("length");

  G4VSolid *sol_module_envelope = new G4Trd("sol_module_envelope",
                                            sin(M_PI / 12.) * rMin, sin(M_PI / 12.) * (rMin + det_height),
                                            segmentlength / 2, segmentlength / 2,
                                            det_height / 2);

  G4LogicalVolume *log_module_envelope = new G4LogicalVolume(sol_module_envelope, Air, "log_module_envelope");

  G4double cooling_plate_height = 2.5 * mm; //6.35
  G4VSolid *sol_cooling_plate = new G4Trd("sol_cooling_plate",
                                          sin(M_PI / 12.) * (rMin + det_height / 2 - cooling_plate_height / 2),
                                          sin(M_PI / 12.) * (rMin + det_height / 2 + cooling_plate_height / 2),
                                          segmentlength / 2, segmentlength / 2,
                                          cooling_plate_height / 2);
  bool doCooling = false;
  if(doCooling){
    G4double diameter_coolingtube = cooling_plate_height-0.5*mm;
    G4double wallthickness_coolingtube = diameter_coolingtube/5;
    G4VSolid *sol_cutout_tube = new G4Cons("sol_cutout_tube",
                                          0, 1.01 * diameter_coolingtube / 2,
                                          0, 1.01 * diameter_coolingtube / 2,
                                          (segmentlength * 1.1) / 2,
                                          0, 2 * M_PI);
    G4VSolid *sol_cooling_tube = new G4Cons("sol_cooling_tube",
                                            (diameter_coolingtube - 2 * wallthickness_coolingtube) / 2, diameter_coolingtube / 2,
                                            (diameter_coolingtube - 2 * wallthickness_coolingtube) / 2, diameter_coolingtube / 2,
                                            (segmentlength - 0.2 * mm) / 2,
                                            0, 2 * M_PI);
    G4LogicalVolume *Log_cooling_tube = new G4LogicalVolume(sol_cooling_tube,  //
                                                            G4Material::GetMaterial("G4_Al"), "Log_cooling_tube");
    RegisterLogicalVolume(Log_cooling_tube);
    m_DisplayAction->AddVolume(Log_cooling_tube, "Cooling_tube");

    G4VSolid *sol_water_cooling = new G4Cons("sol_water_cooling",
                                            0, (diameter_coolingtube - 2 * wallthickness_coolingtube) / 2,
                                            0, (diameter_coolingtube - 2 * wallthickness_coolingtube) / 2,
                                            (segmentlength - 0.3 * mm) / 2,
                                            0, 2 * M_PI);
    G4LogicalVolume *Log_water_cooling = new G4LogicalVolume(sol_water_cooling,  //
                                                            G4Material::GetMaterial("G4_WATER"), "Log_water_cooling");
    RegisterLogicalVolume(Log_water_cooling);
    m_DisplayAction->AddVolume(Log_water_cooling, "Water_cooling");

    G4RotationMatrix *rotcooling = new G4RotationMatrix();
    rotcooling->rotateX(M_PI / 2);
    G4double leftedgeCU = sin(M_PI / 12.) * (rMin + det_height / 2 + cooling_plate_height / 2);
    int maxicup = 12;
    if (rMin < 85 * cm) maxicup = 11;
    if (rMin < 66 * cm) maxicup = 9;
    if (rMin < 55 * cm) maxicup = 7;
    for (int icup = 0; icup < maxicup; icup++)
    {
      G4double edgeshift = 0;
      if (icup == 0) edgeshift = baseplate_width / 4;
      if (icup == (maxicup - 1)) edgeshift = -baseplate_width / 4;
      sol_cooling_plate = new G4SubtractionSolid(G4String("sol_cooling_plate_cu1"), sol_cooling_plate, sol_cutout_tube, rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift, 0, cooling_plate_height / 2 - diameter_coolingtube / 2));
      RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift, 0, cooling_plate_height / 2 - diameter_coolingtube / 2), Log_cooling_tube,
                                              "cooling_tube_Physical_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),
                            false);
      RegisterPhysicalVolume(new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU + icup * 1.5 * baseplate_width + edgeshift, 0, cooling_plate_height / 2 - diameter_coolingtube / 2), Log_water_cooling,
                                              "water_cooling_Physical_" + std::to_string(icup), log_module_envelope, false, 0, overlapcheck_sector),
                            false);
    }
  }
  G4LogicalVolume *log_cooling_plate = new G4LogicalVolume(sol_cooling_plate, G4Material::GetMaterial("G4_Al"), "log_cooling_plate_barrel");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, 0), log_cooling_plate,
                                           "physical_cooling_plate", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterLogicalVolume(log_cooling_plate);
  m_DisplayAction->AddVolume(log_cooling_plate, "CoolingPlate");

  G4double cooling_plate_epoxy_height = 0.08 * mm;
  G4VSolid *sol_cooling_plate_epoxy = new G4Trd("sol_cooling_plate_epoxy",
                                                sin(M_PI / 12.) * (rMin + det_height / 2 + cooling_plate_height / 2),
                                                sin(M_PI / 12.) * (rMin + det_height / 2 + cooling_plate_height / 2 + cooling_plate_epoxy_height),
                                                segmentlength / 2, segmentlength / 2,
                                                cooling_plate_epoxy_height / 2);
  G4LogicalVolume *log_cooling_plate_epoxy = new G4LogicalVolume(sol_cooling_plate_epoxy, mat_Epoxy, "log_cooling_plate_barrel_epoxy");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, cooling_plate_height / 2 + cooling_plate_epoxy_height / 2), log_cooling_plate_epoxy,
                                           "physical_cooling_plate_epoxy", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterLogicalVolume(log_cooling_plate_epoxy);
  m_DisplayAction->AddVolume(log_cooling_plate_epoxy, "Epoxy");

  G4double cooling_plate_cover_height = 0.81 * mm;
  G4VSolid *sol_cooling_plate_cover = new G4Trd("sol_cooling_plate_cover",
                                                sin(M_PI / 12.) * (rMin + det_height / 2 + cooling_plate_height / 2 + cooling_plate_epoxy_height),
                                                sin(M_PI / 12.) * (rMin + det_height / 2 + cooling_plate_height / 2 + cooling_plate_epoxy_height + cooling_plate_cover_height),
                                                segmentlength / 2, segmentlength / 2,
                                                cooling_plate_cover_height / 2);
  G4LogicalVolume *log_cooling_plate_cover = new G4LogicalVolume(sol_cooling_plate_cover, G4Material::GetMaterial("G4_Al"), "log_cooling_plate_barrel_cover");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, cooling_plate_height / 2 + cooling_plate_cover_height / 2 + cooling_plate_epoxy_height), log_cooling_plate_cover,
                                           "physical_cooling_plate_cover", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterLogicalVolume(log_cooling_plate_cover);
  m_DisplayAction->AddVolume(log_cooling_plate_cover, "CoolingPlate");

  // Sensor Module:
  G4double sensor_width = 21.2 * mm;
  G4double sensor_length = 42.0 * mm;
  G4double baseSH_width = baseplate_width / 2;  //-0.15*mm;

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
  G4Material *materialLayer[nLayers] = {G4Material::GetMaterial("G4_GRAPHITE"),
                                        mat_ALN,
                                        G4Material::GetMaterial("G4_GRAPHITE"),
                                        G4Material::GetMaterial("G4_PLEXIGLASS"),
                                        mat_Solder_Tin,
                                        G4Material::GetMaterial("G4_Si"),
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
  G4double offsety = cooling_plate_height / 2 + cooling_plate_epoxy_height + cooling_plate_cover_height + thicknessDet / 2;
  G4double leftedge = sin(M_PI / 12.) * (rMin + det_height / 2 + cooling_plate_height / 2 + cooling_plate_epoxy_height + cooling_plate_cover_height);

  G4RotationMatrix *rotationSensor = new G4RotationMatrix();
  rotationSensor->rotateZ(-M_PI);
  // top side
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t1", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 2 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t2", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 3 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t3", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 5 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t4", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 6 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t5", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 8 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                           "physical_sensor_ladder_t6", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  if (rMin > 55 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 9 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                             "physical_sensor_ladder_t7", log_module_envelope, false, 0, overlapcheck_sector),
                           false);

    RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 11 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                             "physical_sensor_ladder_t8", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  if (rMin > 66 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 12 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                             "physical_sensor_ladder_t9", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  if (rMin > 85 * cm)
  {

    RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 14 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                             "physical_sensor_ladder_t10", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
    RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 15 * baseplate_width + baseplate_width / 2, 0, offsety), log_sensor_ladder,
                                             "physical_sensor_ladder_t11", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  // bottom side
  G4double offsetyDown = cooling_plate_height / 2 + thicknessDet / 2;
  G4RotationMatrix *rotationSensorDown = new G4RotationMatrix();
  rotationSensorDown->rotateY(-M_PI);
  G4RotationMatrix *rotationSensorFlip = new G4RotationMatrix();
  rotationSensorFlip->rotateY(-M_PI);
  rotationSensorFlip->rotateZ(-M_PI);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 1 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b1", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 2 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b2", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 4 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b3", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 5 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b4", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 7 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b5", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 8 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                           "physical_sensor_ladder_b6", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  if (rMin > 55 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 10 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                             "physical_sensor_ladder_b7", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 11 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                             "physical_sensor_ladder_b8", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  if (rMin > 66 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 13 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                             "physical_sensor_ladder_b9", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 14 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
                                             "physical_sensor_ladder_b10", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  if (rMin > 85 * cm)
  {
    // RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 16 * baseplate_width, 0, -offsetyDown), log_sensor_ladder,
    //                                          "physical_sensor_ladder_b11", log_module_envelope, false, 0, overlapcheck_sector),
    //                        false);
  }

  // SERVICE HYBRID
  const int nLayers_SH = 4;
  std::string strLayerName_SH[nLayers_SH] = {
      "ThermalPad",
      "HighSpeedBoard",
      "ConnectorSpace",
      "Powerboard"};
  G4Material *materialLayer_SH[nLayers_SH] = {
      G4Material::GetMaterial("G4_GRAPHITE"),
      G4Material::GetMaterial("G4_POLYSTYRENE"),
      G4Material::GetMaterial("G4_GRAPHITE"),
      G4Material::GetMaterial("G4_POLYSTYRENE")};
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

    G4LogicalVolume *Log_Layer = new G4LogicalVolume(sol_Module_Layer_Raw,  //
                                                     materialLayer_SH[ilay], layer_name + "_Log");
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

  G4double offsety_SH = cooling_plate_height / 2 + cooling_plate_epoxy_height + cooling_plate_cover_height + thicknessDet_SH / 2;
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 1 * baseplate_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t1", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 1 * baseplate_width + 1 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t2", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 3 * baseplate_width + 2 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t3", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 3 * baseplate_width + 3 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t4", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 5 * baseplate_width + 4 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t5", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 5 * baseplate_width + 5 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                           "physical_SH_ladder_t6", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  if (rMin > 55 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 7 * baseplate_width + 6 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                             "physical_SH_ladder_t7", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
    RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 7 * baseplate_width + 7 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                             "physical_SH_ladder_t8", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  if (rMin > 66 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 9 * baseplate_width + 8 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                             "physical_SH_ladder_t9", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
    RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(-leftedge + 9 * baseplate_width + 9 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                             "physical_SH_ladder_t10", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  if (rMin > 85 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensor, G4ThreeVector(-leftedge + 11 * baseplate_width + 10 * baseSH_width + baseSH_width / 2, 0, offsety_SH), log_SH_ladder,
                                             "physical_SH_ladder_t11", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  G4double offsetyDown_SH = cooling_plate_height / 2 + thicknessDet_SH / 2;
  // RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 0 * baseplate_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
  //                                          "physical_SH_ladder_b1", log_module_envelope, false, 0, overlapcheck_sector),
  //                        false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 2 * baseplate_width + 1 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b2", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 2 * baseplate_width + 2 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b3", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 4 * baseplate_width + 3 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b4", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 4 * baseplate_width + 4 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b5", log_module_envelope, false, 0, overlapcheck_sector),
                         false);

  RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 6 * baseplate_width + 5 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                           "physical_SH_ladder_b6", log_module_envelope, false, 0, overlapcheck_sector),
                         false);
  if (rMin > 55 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 6 * baseplate_width + 6 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                             "physical_SH_ladder_b7", log_module_envelope, false, 0, overlapcheck_sector),
                           false);

    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 8 * baseplate_width + 7 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                             "physical_SH_ladder_b8", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  if (rMin > 66 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 8 * baseplate_width + 8 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                             "physical_SH_ladder_b9", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }
  if (rMin > 85 * cm)
  {
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorDown, G4ThreeVector(-leftedge + 10 * baseplate_width + 9 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                             "physical_SH_ladder_b10", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
    RegisterPhysicalVolume(new G4PVPlacement(rotationSensorFlip, G4ThreeVector(-leftedge + 10 * baseplate_width + 10 * baseSH_width + baseSH_width / 2, 0, -offsetyDown_SH), log_SH_ladder,
                                             "physical_SH_ladder_b11", log_module_envelope, false, 0, overlapcheck_sector),
                           false);
  }

  bool doSupport = false;
  G4double support_height = 7 * cm;
  // G4Material* mat_carbonfiber = new G4Material("CarbonFiberSupport", 1.44 * g / cm3, 1);
  // mat_carbonfiber->AddElement(G4Element::GetElement("C"), 1);
  // G4double density;  //z=mean number of protons;
  // G4int ncomponents;
  // carbon+epoxy material
  // G4Material *cfrp_intt = new G4Material("CFRP_INTT", density = 1.69 * g / cm3, ncomponents = 3);
  // cfrp_intt->AddElement(G4Element::GetElement("C"), 10);
  // cfrp_intt->AddElement(G4Element::GetElement("H"), 6);
  // cfrp_intt->AddElement(G4Element::GetElement("O"), 1);

  // SUPPORT STRUCTURES
  G4double support_width = 1 * mm;
  // build components of single segment here
  G4VSolid *Sol_End_Support = new G4Trd("Sol_End_Support",
                                        sin(M_PI / 12.) * (rMin - support_height * 0.9) - 2 * mm, sin(M_PI / 12.) * (rMin) -4 * mm,  // x1, x2
                                        support_width / 2, support_width / 2,                                                        // length
                                        support_height * 0.73 / 2);                                                                  // height

  G4LogicalVolume *Log_End_Support = new G4LogicalVolume(Sol_End_Support, G4Material::GetMaterial("G4_Fe"), "Log_End_Support_Raw");

  // place End side and back side support structure for the segment
  // RegisterPhysicalVolume( new G4PVPlacement(0, G4ThreeVector(0, segmentlength/2-support_width/2, -support_height/2), Log_End_Support,
  //                     "Front_Support_Physical", log_module_envelope, false, 0, overlapcheck_sector), false);
  // RegisterPhysicalVolume( new G4PVPlacement(0, G4ThreeVector(0, -segmentlength/2+support_width/2, -support_height/2), Log_End_Support,
  //                     "Back_Support_Physical", log_module_envelope, false, 0, overlapcheck_sector), false);

  m_DisplayAction->AddVolume(Log_End_Support, "Support");

  // place longitudinal supports left, middle and right side of sector
  G4VSolid *Sol_Longitudinal_Support = new G4Trd("Sol_Longitudinal_Support",
                                                support_width / 2, support_width / 2,                    // x1, x2
                                                segmentlength / 2 - 1 * mm, segmentlength / 2 - 1 * mm,  // length
                                                support_height * 0.73 / 2);                              // height

  G4LogicalVolume *Log_Longitudinal_Support = new G4LogicalVolume(Sol_Longitudinal_Support, G4Material::GetMaterial("G4_Fe"), "Log_Longitudinal_Support_Raw");

  // RegisterPhysicalVolume( new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Log_Longitudinal_Support,
  //                     "Mother_Segment_Raw_Physical_Center", log_module_envelope, false, 0, overlapcheck_sector), false);

  if(doSupport){
    G4RotationMatrix *supportrot = new G4RotationMatrix();
    supportrot->rotateY(-M_PI / 12.);
    if (rMin < 85 * cm)
    {
      RegisterPhysicalVolume(new G4PVPlacement(supportrot, G4ThreeVector(sin(M_PI / 12.) * (rMin - support_height / 2) - support_width / 2, 0, -support_height / 2), Log_Longitudinal_Support,
                                              "Mother_Segment_Raw_Physical_Left", log_module_envelope, false, 0, overlapcheck_sector),
                            false);
      G4RotationMatrix *supportrot2 = new G4RotationMatrix();
      supportrot2->rotateY(M_PI / 12.);
      RegisterPhysicalVolume(new G4PVPlacement(supportrot2, G4ThreeVector(-sin(M_PI / 12.) * (rMin - support_height / 2) + support_width / 2, 0, -support_height / 2), Log_Longitudinal_Support,
                                              "Mother_Segment_Raw_Physical_Right", log_module_envelope, false, 0, overlapcheck_sector),
                            false);
    }
    m_DisplayAction->AddVolume(Log_Longitudinal_Support, "Support");
  }
  RegisterLogicalVolume(log_module_envelope);
  m_DisplayAction->AddVolume(log_module_envelope, "ModuleEnvelope");
  G4double modulesep = 1 * mm;
  G4double moduleShift = -8 * mm;
  if (rMin < 85 * cm) moduleShift = -3 * mm;
  if (rMin < 66 * cm) moduleShift = -1 * mm;
  if (rMin < 55 * cm) moduleShift = 4 * mm;

  for (int isec = 0; isec < 12; isec++)
  {
    // if(isec!=3 && isec!=4)continue; // NOTE REMOVE
    // if(isec!=3)continue; // NOTE REMOVE
    G4RotationMatrix *motherrot = new G4RotationMatrix();
    motherrot->rotateX(M_PI / 2);
    motherrot->rotateY((isec - 3) * 2 * M_PI / 12.);
    // // central segments
    RegisterPhysicalVolume(new G4PVPlacement(motherrot, G4ThreeVector((rMin - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), 0 * modulesep), log_module_envelope,
                                             "Mother_Segment_Raw_Physical_Center_" + std::to_string(isec), DetectorLog_Det, false, 0, overlapcheck_sector),
                           false);
    for (int ilen = 1; ilen < ((detlength / 2 - segmentlength / 2) / segmentlength); ilen++)
    {
      if(doSupport){
        G4RotationMatrix *supfinalrot = new G4RotationMatrix();
        // supfinalrot->rotateX(M_PI/2);
        supfinalrot->rotateX(M_PI / 2);
        supfinalrot->rotateY((isec - 3) * 2 * M_PI / 12.);
        if (ilen == 2 || (ilen == 7))
        {
          if (rMin < 85 * cm)
          {
            RegisterPhysicalVolume(new G4PVPlacement(supfinalrot, G4ThreeVector((rMin - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), ilen * segmentlength + segmentlength / 2 + ilen * modulesep), Log_End_Support,
                                                    "Front_Support_Physical_1_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector),
                                  false);
            RegisterPhysicalVolume(new G4PVPlacement(supfinalrot, G4ThreeVector((rMin - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), -(ilen * segmentlength + segmentlength / 2 + ilen * modulesep)), Log_End_Support,
                                                    "Front_Support_Physical_2_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector),
                                  false);
          }
        }
      }
      // forward segments
      RegisterPhysicalVolume(new G4PVPlacement(motherrot, G4ThreeVector((rMin - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), ilen * segmentlength + ilen * modulesep), log_module_envelope,
                                               "Mother_Segment_Raw_Physical_Fwd_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector),
                             false);
      // backward segments
      RegisterPhysicalVolume(new G4PVPlacement(motherrot, G4ThreeVector((rMin - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), -ilen * segmentlength - ilen * modulesep), log_module_envelope,
                                               "Mother_Segment_Raw_Physical_Bwd_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector),
                             false);
    }
  }
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
  G4Material *mat_Epoxy = G4Material::GetMaterial("EpoxyTTL");
  if (!mat_Epoxy)
  {
    mat_Epoxy = new G4Material("EpoxyTTL", density = 1.16 * g / cm3, natoms = 4);
    mat_Epoxy->AddElement(elH, 32);  // Hydrogen
    mat_Epoxy->AddElement(elN, 2);   // Nitrogen
    mat_Epoxy->AddElement(elO, 4);   // Oxygen
    mat_Epoxy->AddElement(elC, 15);  // Carbon
    // G4Material *mat_Epoxy = G4Material::GetMaterial("Epoxy");
  }
  G4Material *mat_ALN = G4Material::GetMaterial("AluminiumNitrate");
  if (!mat_ALN)
  {
    mat_ALN = new G4Material("AluminiumNitrate", density = 3.255 * g / cm3, ncomponents = 2);
    // G4Material *mat_ALN = new G4Material("AluminiumNitrate", density = 3.255 * g / cm3, ncomponents = 2);
    mat_ALN->AddElement(G4Element::GetElement("Al"), 1);
    mat_ALN->AddElement(G4Element::GetElement("N"), 1);
  }
  G4Material *mat_Solder_Tin = G4Material::GetMaterial("Tin");
  if (!mat_Solder_Tin)
  {
    mat_Solder_Tin = new G4Material("Tin", z = 50., a = 118.7 * g / mole, density = 7.310 * g / cm3);
  }

  G4double det_height = 1.8 * cm;

  //Create the envelope = 'world volume' for the calorimeter
  G4Material *Air = G4Material::GetMaterial("G4_AIR");

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

  G4double cooling_plate_height = 3.35 * mm;
  G4double cooling_plate_epoxy_height = 0.08 * mm;
  G4double cooling_plate_cover_height = 0.81 * mm;

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
  // cooling lines still need to be added TODO
  if (0)
  {
    // TODO
    // G4double diameter_coolingtube = 5*mm;
    // G4double wallthickness_coolingtube = 1 * mm;
    // G4VSolid* sol_cutout_tube = new G4Cons("sol_cutout_tube",
    //                             0, 1.01*diameter_coolingtube / 2,
    //                             0, 1.01*diameter_coolingtube / 2,
    //                             (segmentlength*1.1)/2,
    //                             0, 2 * M_PI);
    // G4VSolid* sol_cooling_tube = new G4Cons("sol_cooling_tube",
    //                             (diameter_coolingtube - 2*wallthickness_coolingtube) / 2, diameter_coolingtube / 2,
    //                             (diameter_coolingtube - 2*wallthickness_coolingtube) / 2, diameter_coolingtube / 2,
    //                             (segmentlength-0.2*mm)/2,
    //                             0, 2 * M_PI);
    // G4LogicalVolume *Log_cooling_tube = new G4LogicalVolume(sol_cooling_tube,  //
    //                                                       G4Material::GetMaterial("G4_Al"), "Log_cooling_tube");
    // RegisterLogicalVolume(Log_cooling_tube);
    // m_DisplayAction->AddVolume(Log_cooling_tube, "Cooling_tube");

    // G4VSolid* sol_water_cooling = new G4Cons("sol_water_cooling",
    //                             0, (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
    //                             0, (diameter_coolingtube - 2*wallthickness_coolingtube) / 2,
    //                             (segmentlength-0.3*mm)/2,
    //                             0, 2 * M_PI);
    // G4LogicalVolume *Log_water_cooling = new G4LogicalVolume(sol_water_cooling,  //
    //                                                       G4Material::GetMaterial("G4_WATER"), "Log_water_cooling");
    // RegisterLogicalVolume(Log_water_cooling);
    // m_DisplayAction->AddVolume(Log_water_cooling, "Water_cooling");

    // G4RotationMatrix *rotcooling = new G4RotationMatrix();
    // rotcooling->rotateX(M_PI/2);
    // G4double leftedgeCU = sin(M_PI/12.)*(rMin+det_height/2+cooling_plate_height/2);
    // for(int icup=0;icup<11;icup++){
    //   G4double edgeshift = 0;
    //   if(icup==0) edgeshift = baseplate_width/4;
    //   if(icup==10) edgeshift = -baseplate_width/4;
    //   sol_cooling_plate = new G4SubtractionSolid(G4String("sol_cooling_plate_cu1"), sol_cooling_plate, sol_cutout_tube
    //                                                     ,rotcooling ,G4ThreeVector( -leftedgeCU+icup*1.5*baseplate_width+edgeshift , 0 ,cooling_plate_height/2 - diameter_coolingtube/2));
    //       RegisterPhysicalVolume( new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU+icup*1.5*baseplate_width+edgeshift, 0, cooling_plate_height/2 - diameter_coolingtube/2 ), Log_cooling_tube,
    //                             "cooling_tube_Physical_" + icup, log_module_envelope, false, 0, overlapcheck_sector), false);
    //       RegisterPhysicalVolume( new G4PVPlacement(rotcooling, G4ThreeVector(-leftedgeCU+icup*1.5*baseplate_width+edgeshift, 0, cooling_plate_height/2 - diameter_coolingtube/2 ), Log_water_cooling,
    //                             "water_cooling_Physical_" + icup, log_module_envelope, false, 0, overlapcheck_sector), false);
    // }
  }

  G4LogicalVolume *log_cooling_plate = new G4LogicalVolume(sol_cooling_plate, G4Material::GetMaterial("G4_Al"), "log_cooling_plate_fwd");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, 0), log_cooling_plate,
                                            "physical_cooling_plate", log_module_envelope, false, 0, overlapcheck_sector),
                          false);
  RegisterLogicalVolume(log_cooling_plate);
  m_DisplayAction->AddVolume(log_cooling_plate, "CoolingPlate");

  G4VSolid *sol_cooling_plate_epoxy = new G4Cons("ttlmodule_beampipe_cutout_cooling_plate_epoxy",
                                                  0, rMax,
                                                  0, rMax,
                                                  cooling_plate_epoxy_height / 2.0,
                                                  0, 2 * M_PI);
  sol_cooling_plate_epoxy = new G4SubtractionSolid(G4String("sol_cooling_plate_epoxy"), sol_cooling_plate_epoxy, beampipe_cutout, 0, G4ThreeVector(xoffset, 0, 0.));
  G4LogicalVolume *log_cooling_plate_epoxy = new G4LogicalVolume(sol_cooling_plate_epoxy, mat_Epoxy, "log_cooling_plate_fwd_epoxy");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, cooling_plate_height / 2 + cooling_plate_epoxy_height / 2), log_cooling_plate_epoxy,
                                            "physical_cooling_plate_epoxy", log_module_envelope, false, 0, overlapcheck_sector),
                          false);
  RegisterLogicalVolume(log_cooling_plate_epoxy);
  m_DisplayAction->AddVolume(log_cooling_plate_epoxy, "Epoxy");

  G4VSolid *sol_cooling_plate_cover = new G4Cons("ttlmodule_beampipe_sol_cooling_plate_cover",
                                                  0, rMax,
                                                  0, rMax,
                                                  cooling_plate_cover_height / 2.0,
                                                  0, 2 * M_PI);
  sol_cooling_plate_cover = new G4SubtractionSolid(G4String("sol_cooling_plate_cover"), sol_cooling_plate_cover, beampipe_cutout, 0, G4ThreeVector(xoffset, 0, 0.));
  G4LogicalVolume *log_cooling_plate_cover = new G4LogicalVolume(sol_cooling_plate_cover, G4Material::GetMaterial("G4_Al"), "log_cooling_plate_fwd_cover");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, cooling_plate_height / 2 + cooling_plate_cover_height / 2 + cooling_plate_epoxy_height), log_cooling_plate_cover,
                                            "physical_cooling_plate_cover", log_module_envelope, false, 0, overlapcheck_sector),
                          false);
  RegisterLogicalVolume(log_cooling_plate_cover);
  m_DisplayAction->AddVolume(log_cooling_plate_cover, "CoolingPlate");



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
  G4Material *materialLayer[nLayers] = {G4Material::GetMaterial("G4_GRAPHITE"), mat_ALN, G4Material::GetMaterial("G4_GRAPHITE"), G4Material::GetMaterial("G4_PLEXIGLASS"), mat_Solder_Tin, G4Material::GetMaterial("G4_Si"), mat_Epoxy, mat_ALN};
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
      G4Material::GetMaterial("G4_GRAPHITE"),
      G4Material::GetMaterial("G4_POLYSTYRENE"),
      G4Material::GetMaterial("G4_GRAPHITE"),
      G4Material::GetMaterial("G4_POLYSTYRENE")};
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

  G4double offsetzFront = cooling_plate_height / 2 + cooling_plate_epoxy_height + cooling_plate_cover_height + thicknessDet_SH / 2;
  G4double offsetzBack = -cooling_plate_height / 2 - thicknessDet_SH / 2;
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

        // create mother volume with space for numSensorsRow towers along x-axis
        auto TTLDetRowRightSolid    = new G4Box("TTLDetRowRightBox" + std::to_string(row), (((numSensorsRow-1) /2 - numSensorRightAdd)) * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        auto TTLDetRowRightLogical  = new G4LogicalVolume(TTLDetRowRightSolid,Air,"TTLDetRowRightLogical" + std::to_string(row));
        // replicate singletower tower design numSensorsRow times along x-axis
        new G4PVReplica("TTLDetRowRightPhysical" + std::to_string(row),log_sensor_and_readout,TTLDetRowRightLogical,
                        kXAxis,((numSensorsRow-1) /2 - numSensorRightAdd ),segmentlength);

        new G4PVPlacement(row%2==0 ? rotationSensor : 0, G4ThreeVector(( (((numSensorsRow-1) /2 - numSensorRightAdd)) * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), (row*fullsensor_width), offsetzFront),
                      TTLDetRowRightLogical, "TTLDetRowRightPlacedFront" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector);

        new G4PVPlacement(row%2==0 ? rotationSensorBck2 : rotationSensorBck1, G4ThreeVector(( (((numSensorsRow-1) /2 - numSensorRightAdd)) * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), (row*fullsensor_width), offsetzBack),
                      TTLDetRowRightLogical, "TTLDetRowRightPlacedBack" + std::to_string(row), log_module_envelope, 0, false, overlapcheck_sector);


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
      }
    } else {
        // create mother volume with space for numSensorsRow towers along x-axis
        auto TTLDetRowSolid    = new G4Box("TTLDetRowBox" + std::to_string(row), numSensorsRow * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        auto TTLDetRowLogical  = new G4LogicalVolume(TTLDetRowSolid,Air,"TTLDetRowLogical" + std::to_string(row));
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

        // for(int icol=-(numSensorsRow-1)/2;icol<=(numSensorsRow-1)/2;icol++){
        //   RegisterPhysicalVolume(new G4PVPlacement(row%2==0 ? rotationSensor : 0, G4ThreeVector(icol*segmentlength, (row*fullsensor_width), offsetzFront),
        //                 log_sensor_and_readout, "TTLDetRowPhysicalFront_j" + std::to_string(row) + "_k"  + std::to_string(icol), log_module_envelope, false, 0, overlapcheck_sector),false);

        //   RegisterPhysicalVolume(new G4PVPlacement(row%2==0 ? rotationSensorBck2 : rotationSensorBck1, G4ThreeVector(icol*segmentlength, (row*fullsensor_width), offsetzBack),
        //                 log_sensor_and_readout, "TTLDetRowPhysicalBack_j" + std::to_string(row) + "_k"  + std::to_string(icol), log_module_envelope, false, 0, overlapcheck_sector),false);
        // }
    }
  }
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
