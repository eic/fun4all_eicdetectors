#include "PHG4HybridHomogeneousCalorimeterDetector.h"
#include "PHG4HybridHomogeneousCalorimeterDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Element.hh>  // for G4Element
#include <Geant4/G4LogicalBorderSurface.hh>
#include <Geant4/G4LogicalSkinSurface.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MaterialPropertiesTable.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Polyhedra.hh>       // for G4RotationMatrix
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Types.hh>        // for G4double
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________________
PHG4HybridHomogeneousCalorimeterDetector::PHG4HybridHomogeneousCalorimeterDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_SuperDetector("NONE")
  , m_Params(parameters)
  , m_DisplayAction(dynamic_cast<PHG4HybridHomogeneousCalorimeterDisplayAction*>(subsys->GetDisplayAction()))
  , _towerlogicnameprefix("HybridHomogeneousCalorimeterTower")
  , m_IsActive(m_Params->get_int_param("active"))
  , m_AbsorberActive(m_Params->get_int_param("absorberactive"))
  , m_doLightProp(false)
{
}

//_______________________________________________________________________
int PHG4HybridHomogeneousCalorimeterDetector::IsInCrystalCalorimeter(G4VPhysicalVolume* volume) const
{
  if (m_IsActive)
  {
    if (m_ActiveVolumeSet.find(volume) != m_ActiveVolumeSet.end())
    {
      return GetCaloType();
    }
  }
  if (m_AbsorberActive)
  {
    if (m_PassiveVolumeSet.find(volume) != m_PassiveVolumeSet.end())
    {
      return -1;
    }
  }
  return 0;
}

//_______________________________________________________________________
void PHG4HybridHomogeneousCalorimeterDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4HybridHomogeneousCalorimeterDetector: Begin Construction" << endl;
  }

  if (m_Params->get_string_param("mappingtower").empty())
  {
    cout << "ERROR in PHG4HybridHomogeneousCalorimeterDetector: No tower mapping file specified. Abort detector construction." << endl;
    cout << "Please run set_string_param(\"mappingtower\", std::string filename ) first." << endl;
    exit(1);
  }

  /* Read parameters for detector construction and mapping from file */
  ParseParametersFromTable();

  // G4LogicalVolume* support_frame = ConstructSupportFrame(logicWorld);

  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  /* Place calorimeter tower within envelope */
  //  PlaceTower(eemc_envelope_log, singletower);
  PlaceTower(logicWorld, singletower);  //,support_frame);

  return;
}
//_______________________________________________________________________
G4LogicalVolume* PHG4HybridHomogeneousCalorimeterDetector::ConstructSupportFrame(G4LogicalVolume* eemcenvelope)
{
  G4double frame_r_in = m_Params->get_double_param("rMin2") * cm;
  G4double frame_r_out = m_Params->get_double_param("rMax2") * cm;
  G4double frame_dz = m_Params->get_double_param("dz") * cm;
  G4double frame_z_pos = m_Params->get_double_param("place_z") * cm;

  G4int nSides = 12;
  const G4int nZSlices = 4;
  G4double zPosSlice[nZSlices] = {0, 0.001, frame_dz - 0.001, frame_dz};
  G4double rinPosSlice[nZSlices] = {frame_r_in, frame_r_in, frame_r_in, frame_r_in};
  G4double routPosSlice[nZSlices] = {frame_r_out, frame_r_out, frame_r_out, frame_r_out};

  G4VSolid* frame_full_solid = new G4Polyhedra(G4String("frame_full_solid"),
                                               0,
                                               2 * M_PI * rad,
                                               nSides,
                                               nZSlices,
                                               zPosSlice,
                                               rinPosSlice,
                                               routPosSlice);
  G4LogicalVolume* frame_full_logic = new G4LogicalVolume(frame_full_solid, G4Material::GetMaterial("G4_Fe"), "frame_full_logic", 0, 0, 0);
  GetDisplayAction()->AddVolume(frame_full_logic, "WIP");

  G4RotationMatrix* rotFrame = new G4RotationMatrix();
  rotFrame->rotateZ(M_PI / 12);
  new G4PVPlacement(rotFrame, G4ThreeVector(0, 0, frame_z_pos - frame_dz / 2),
                    frame_full_logic,
                    "frame_full_placed",
                    eemcenvelope,
                    0, 0, OverlapCheck());

  return frame_full_logic;
}

//_______________________________________________________________________
G4LogicalVolume* PHG4HybridHomogeneousCalorimeterDetector::ConstructTower()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4HybridHomogeneousCalorimeterDetector: Build logical volume for single tower..." << endl;
  }

  /* dimensions of full tower */
  G4double crystal_dx = m_Params->get_double_param("crystal_dx") * cm;
  G4double crystal_dy = m_Params->get_double_param("crystal_dy") * cm;
  G4double crystal_dz = m_Params->get_double_param("crystal_dz") * cm;

  G4double carbon_thickness = m_Params->get_double_param("carbon_gap") * cm;
  G4double airgap_crystal_carbon = m_Params->get_double_param("air_gap") * cm;

  G4double reflective_foil_thickness = m_Params->get_double_param("reflective_foil_thickness") * cm;
  G4double tedlar_thickness = m_Params->get_double_param("tedlar_thickness") * cm;
  bool doWrapping = false;
  if (reflective_foil_thickness > 0 || tedlar_thickness > 0) doWrapping = true;

  G4int sensor_count = m_Params->get_int_param("sensor_count");
  G4double sensor_dimension = m_Params->get_double_param("sensor_dimension") * cm;
  G4double sensor_thickness = m_Params->get_double_param("sensor_thickness") * cm;
  bool doSensors = false;
  if (sensor_dimension > 0 && sensor_count > 0 && sensor_thickness > 0) doSensors = true;

  G4double tower_dx = crystal_dx + 2 * (carbon_thickness + airgap_crystal_carbon + reflective_foil_thickness + tedlar_thickness);
  G4double tower_dy = crystal_dy + 2 * (carbon_thickness + airgap_crystal_carbon + reflective_foil_thickness + tedlar_thickness);
  G4double tower_dz = crystal_dz + 2 * (carbon_thickness);
  if (doSensors) tower_dz = crystal_dz + 2 * (carbon_thickness) + sensor_thickness;
  int carbon_frame_style = m_Params->get_int_param("carbon_frame_style");

  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));

  /* create logical volume for single tower */
  // Building the single tower mother volume first
  // Then the crystal/sci-glass will be put in this volume
  // The shell will also placed in. The rest space leave for the air gap

  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"), tower_dx / 2.0, tower_dy / 2.0, tower_dz / 2.0);
  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid, WorldMaterial, "single_tower_logic", 0, 0, 0);
  GetDisplayAction()->AddVolume(single_tower_logic, "Invisible");
  /* create geometry volume for crystal inside single_tower */
  G4VSolid* solid_crystal = new G4Box(G4String("single_crystal_solid"),
                                      crystal_dx / 2.0,
                                      crystal_dy / 2.0,
                                      crystal_dz / 2.0);

  G4Material* material_shell = GetCarbonFiber();
  if (carbon_frame_style == 0 && carbon_thickness > 0)
  {
    /* create geometry volume for frame (carbon fiber shell) inside single_tower */
    G4VSolid* Carbon_hunk_solid = new G4Box(G4String("Carbon_hunk_solid"),
                                            tower_dx / 2.0,
                                            tower_dy / 2.0,
                                            ((tower_dz / 2.0) - 1 * mm));

    G4VSolid* cutout_solid = new G4Box(G4String("lead_solid"), crystal_dx / 2.0, crystal_dy / 2.0, crystal_dz);

    G4SubtractionSolid* Carbon_shell_solid = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
                                                                    Carbon_hunk_solid,
                                                                    cutout_solid,
                                                                    0,
                                                                    G4ThreeVector(0.00 * mm, 0.00 * mm, 0.00 * mm));

    G4LogicalVolume* Carbon_shell_logic = new G4LogicalVolume(Carbon_shell_solid, material_shell, "Carbon_shell_logic", 0, 0, 0);
    GetDisplayAction()->AddVolume(Carbon_shell_logic, "CarbonShell");

    /* Place structural frame in logical tower volume */
    // ostringstream name_shell;
    // name_shell.str("");
    // name_shell << _towerlogicnameprefix << "_single_absorber";
    string name_shell = _towerlogicnameprefix + "_single_shell";
    G4VPhysicalVolume* physvol_carbon = new G4PVPlacement(0, G4ThreeVector(0, 0, sensor_thickness / 2), Carbon_shell_logic, name_shell, single_tower_logic, 0, 0, OverlapCheck());
    m_PassiveVolumeSet.insert(physvol_carbon);
  }
  else if (carbon_frame_style == 1)
  {
    // the carbon spacers should go a few cm deep between the crystals
    G4double carbon_frame_depth = m_Params->get_double_param("carbon_frame_depth") * cm;
    G4VSolid* Carbon_hunk_solid = new G4Box(G4String("Carbon_hunk_solid"),
                                            tower_dx / 2.0,
                                            tower_dy / 2.0,
                                            (carbon_frame_depth / 2.0));
    G4VSolid* cutout_solid = new G4Box(G4String("cutout_solid"), (tower_dx - 2 * (carbon_thickness)) / 2.0, (tower_dy - 2 * (carbon_thickness)) / 2.0, crystal_dz);

    G4SubtractionSolid* Carbon_gap_solid = new G4SubtractionSolid(G4String("Carbon_gap_solid"),
                                                                  Carbon_hunk_solid,
                                                                  cutout_solid,
                                                                  0,
                                                                  G4ThreeVector(0.00 * mm, 0.00 * mm, 0.00 * mm));
    G4LogicalVolume* logic_gap = new G4LogicalVolume(Carbon_gap_solid, material_shell, "Carbon_gap_logic", 0, 0, 0);
    GetDisplayAction()->AddVolume(logic_gap, "CarbonShell");

    string name_gap = _towerlogicnameprefix + "_single_shell";
    G4VPhysicalVolume* physvol_carbon_gap_front = new G4PVPlacement(0, G4ThreeVector(0, 0, tower_dz / 2 - carbon_frame_depth / 2 - carbon_thickness), logic_gap, name_gap + "_gap_front", single_tower_logic, 0, 0, OverlapCheck());
    m_PassiveVolumeSet.insert(physvol_carbon_gap_front);
    G4VPhysicalVolume* physvol_carbon_gap_back = new G4PVPlacement(0, G4ThreeVector(0, 0, -tower_dz / 2 + carbon_frame_depth / 2 + carbon_thickness + sensor_thickness), logic_gap, name_gap + "_gap_back", single_tower_logic, 0, 0, OverlapCheck());
    m_PassiveVolumeSet.insert(physvol_carbon_gap_back);

    G4double carbon_face_lip = m_Params->get_double_param("carbon_face_lip") * cm;

    G4VSolid* Carbon_hunk_solid_face = new G4Box(G4String("Carbon_hunk_solid_face"),
                                                 tower_dx / 2.0,
                                                 tower_dy / 2.0,
                                                 (carbon_thickness / 2.0));
    G4VSolid* cutout_solid_face = new G4Box(G4String("cutout_solid_face"), (tower_dx - 2.0 * carbon_face_lip) / 2.0, (tower_dy - 2.0 * carbon_face_lip) / 2.0, (tower_dz));

    G4SubtractionSolid* Carbon_face_solid = new G4SubtractionSolid(G4String("Carbon_face_solid"),
                                                                   Carbon_hunk_solid_face,
                                                                   cutout_solid_face,
                                                                   0,
                                                                   G4ThreeVector(0.00 * mm, 0.00 * mm, 0.00 * mm));
    G4LogicalVolume* logic_face = new G4LogicalVolume(Carbon_face_solid, material_shell, "Carbon_face_logic", 0, 0, 0);
    GetDisplayAction()->AddVolume(logic_face, "CarbonShell");

    string name_face = _towerlogicnameprefix + "_single_shell";
    G4VPhysicalVolume* physvol_carbon_face_front = new G4PVPlacement(0, G4ThreeVector(0, 0, tower_dz / 2 - carbon_thickness / 2), logic_face, name_face + "_face_front", single_tower_logic, 0, 0, OverlapCheck());
    m_PassiveVolumeSet.insert(physvol_carbon_face_front);
    G4VPhysicalVolume* physvol_carbon_face_back = new G4PVPlacement(0, G4ThreeVector(0, 0, -tower_dz / 2 + carbon_thickness / 2 + sensor_thickness), logic_face, name_face + "_face_back", single_tower_logic, 0, 0, OverlapCheck());
    m_PassiveVolumeSet.insert(physvol_carbon_face_back);
  }
  // wrap towers in VM2000 and Tedlar foils if requested
  if (doWrapping)
  {
    G4VSolid* VM2000_hunk_solid = new G4Box(G4String("VM2000_hunk_solid"),
                                            (crystal_dx + 2 * reflective_foil_thickness) / 2.0,
                                            (crystal_dy + 2 * reflective_foil_thickness) / 2.0,
                                            ((crystal_dz / 2.0) - 1 * mm));

    G4VSolid* cutout_solid_VM2000 = new G4Box(G4String("cutout_solid_VM2000"), crystal_dx / 2.0, crystal_dy / 2.0, crystal_dz);

    G4SubtractionSolid* VM2000_foil_solid = new G4SubtractionSolid(G4String("VM2000_foil_solid"),
                                                                   VM2000_hunk_solid,
                                                                   cutout_solid_VM2000,
                                                                   0,
                                                                   G4ThreeVector(0.00 * mm, 0.00 * mm, 0.00 * mm));

    G4LogicalVolume* VM2000_foil_logic = new G4LogicalVolume(VM2000_foil_solid, GetVM2000Material(), "VM2000_foil_logic", 0, 0, 0);
    GetDisplayAction()->AddVolume(VM2000_foil_logic, "VM2000");

    string name_VM2000foil = _towerlogicnameprefix + "_single_foil_VM2000";
    G4VPhysicalVolume* physvol_VM2000 = new G4PVPlacement(0, G4ThreeVector(0, 0, sensor_thickness / 2), VM2000_foil_logic, name_VM2000foil, single_tower_logic, 0, 0, OverlapCheck());
    m_PassiveVolumeSet.insert(physvol_VM2000);

    /* create geometry volume for frame (carbon fiber shell) inside single_tower */
    G4VSolid* Tedlar_hunk_solid = new G4Box(G4String("Tedlar_hunk_solid"),
                                            (crystal_dx + 2 * reflective_foil_thickness + 2 * tedlar_thickness) / 2.0,
                                            (crystal_dy + 2 * reflective_foil_thickness + 2 * tedlar_thickness) / 2.0,
                                            ((crystal_dz / 2.0) - 1 * mm));

    G4VSolid* cutout_solid_Tedlar = new G4Box(G4String("cutout_solid_Tedlar"), (crystal_dx + 2 * reflective_foil_thickness) / 2.0, (crystal_dy + 2 * reflective_foil_thickness) / 2.0, crystal_dz);

    G4SubtractionSolid* Tedlar_foil_solid = new G4SubtractionSolid(G4String("Tedlar_foil_solid"),
                                                                   Tedlar_hunk_solid,
                                                                   cutout_solid_Tedlar,
                                                                   0,
                                                                   G4ThreeVector(0.00 * mm, 0.00 * mm, 0.00 * mm));

    G4LogicalVolume* Tedlar_foil_logic = new G4LogicalVolume(Tedlar_foil_solid, GetTedlarMaterial(), "Tedlar_foil_logic", 0, 0, 0);
    GetDisplayAction()->AddVolume(Tedlar_foil_logic, "Tedlar");

    string name_Tedlarfoil = _towerlogicnameprefix + "_single_foil_Tedlar";
    G4VPhysicalVolume* physvol_Tedlar = new G4PVPlacement(0, G4ThreeVector(0, 0, sensor_thickness / 2), Tedlar_foil_logic, name_Tedlarfoil, single_tower_logic, 0, 0, OverlapCheck());
    m_PassiveVolumeSet.insert(physvol_Tedlar);
  }
  /* create logical volumes for crystal inside single_tower */
  G4double M_para = m_Params->get_double_param("material");
  // set default material for hom calo
  G4Material* material_Scin = G4Material::GetMaterial("G4_PbWO4");
  if (M_para > 0) material_Scin = GetScintillatorMaterial(M_para);

  if (m_doLightProp)
  {
    //scintillation and optical properties for the crystal
    CrystalTable(material_Scin);
  }

  G4LogicalVolume* logic_crystal = new G4LogicalVolume(solid_crystal, material_Scin, "single_crystal_logic", 0, 0, 0);

  if (m_doLightProp)
  {
    //crystal optical surface
    SurfaceTable(logic_crystal);
  }
  GetDisplayAction()->AddVolume(logic_crystal, "Crystal");

  /* Place crystal in logical tower volume */
  string name_crystal = _towerlogicnameprefix + "_single_crystal";
  G4VPhysicalVolume* physvol_crys = new G4PVPlacement(0, G4ThreeVector(0, 0, sensor_thickness / 2), logic_crystal, name_crystal, single_tower_logic, 0, 0, OverlapCheck());
  m_ActiveVolumeSet.insert(physvol_crys);

  if (doSensors)
  {
    G4VSolid* single_sensor_solid = new G4Box("single_sensor_solid", sensor_dimension / 2., sensor_dimension / 2., sensor_thickness / 2.);
    G4Material* material_Sensor = G4Material::GetMaterial("G4_Si");

    G4LogicalVolume* single_sensor_logic = new G4LogicalVolume(single_sensor_solid, material_Sensor, "single_sensor_logic", 0, 0, 0);
    GetDisplayAction()->AddVolume(single_sensor_logic, "Sensor");

    string name_sensor = _towerlogicnameprefix + "_single_sensor";

    G4VPhysicalVolume* physvol_sensor_0 = new G4PVPlacement(0, G4ThreeVector(0, 0, -tower_dz / 2 + carbon_thickness + sensor_thickness / 2), single_sensor_logic, name_sensor, single_tower_logic, 0, 0, OverlapCheck());
    m_ActiveVolumeSet.insert(physvol_sensor_0);
    if (m_doLightProp)
    {
      MakeBoundary(physvol_crys, physvol_sensor_0);
    }
  }
  if (Verbosity() > 0)
    cout << "PHG4HybridHomogeneousCalorimeterDetector: Building logical volume for single tower done." << endl;

  return single_tower_logic;
}

int PHG4HybridHomogeneousCalorimeterDetector::PlaceTower(G4LogicalVolume* eemcenvelope, G4LogicalVolume* singletower)  //, G4LogicalVolume* framelogic)
{
  /* Loop over all tower positions in vector and place tower */
  typedef std::map<std::string, towerposition>::iterator it_type;
  for (it_type iterator = _map_tower.begin(); iterator != _map_tower.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4HybridHomogeneousCalorimeterDetector: Place tower " << iterator->first
           << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k
           << " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << endl;
    }
    int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;
    new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                      singletower,
                      iterator->first,
                      eemcenvelope,
                      0, copyno, OverlapCheck());
  }

  return 0;
}

//_______________________________________________________________________
G4Material*
PHG4HybridHomogeneousCalorimeterDetector::GetScintillatorMaterial(float setting)
{
  G4Element* ele_O = new G4Element("Oxygen", "O", 8., 16.00 * g / mole);
  G4Element* ele_Si = new G4Element("Silicon", "Si", 14., 28.09 * g / mole);
  G4Element* ele_B = new G4Element("Boron", "B", 5., 10.811 * g / mole);
  G4Element* ele_Na = new G4Element("Sodium", "Na", 11., 22.99 * g / mole);
  G4Element* ele_Mg = new G4Element("Magnesium", "Mg", 12., 24.30 * g / mole);
  G4Element* ele_Pb = new G4Element("Lead", "Pb", 82., 207.2 * g / mole);
  G4Element* ele_Ba = new G4Element("Barium", "Ba", 56., 137.3 * g / mole);
  G4Element* ele_Gd = new G4Element("Gadolinium", "Gd", 64., 157.3 * g / mole);

  G4Material* material_Scin = G4Material::GetMaterial("G4_PbWO4");

  if ((setting > 0.) && (setting < 1.))
  {
    material_Scin = G4Material::GetMaterial("G4_PbWO4");
    // g4MatData.push_back(0.0333333*mm/MeV);
    if (Verbosity()) cout << "Set G4_PbWO4..." << endl;
  }
  else if ((setting > 1.) && (setting < 2.))
  {
    material_Scin = G4Material::GetMaterial("G4_GLASS_LEAD");
    if (Verbosity()) cout << "Set G4_GLASS_LEAD..." << endl;
  }
  else if ((setting > 2.) && (setting < 3.))
  {
    material_Scin = G4Material::GetMaterial("G4_BARIUM_SULFATE");
    if (Verbosity()) cout << "Set G4_BARIUM_SULFATE..." << endl;
  }
  else if ((setting > 3.) && (setting < 4.))
  {
    material_Scin = G4Material::GetMaterial("G4_CESIUM_IODIDE");
    if (Verbosity()) cout << "Set G4_CESIUM_IODIDE..." << endl;
  }
  else if ((setting > 4.) && (setting < 5.))
  {
    material_Scin = new G4Material("material_Scin", 4.5 * g / cm3, 5);
    material_Scin->AddElement(ele_Si, 21.9 * perCent);
    material_Scin->AddElement(ele_B, 8.8 * perCent);
    material_Scin->AddElement(ele_Na, 10.4 * perCent);
    material_Scin->AddElement(ele_Mg, 6.5 * perCent);
    material_Scin->AddElement(ele_O, 52.4 * perCent);

    if (Verbosity()) cout << "Set Sciglass..." << endl;
  }
  else if ((setting > 5.) && (setting < 6.))
  {
    material_Scin = new G4Material("material_Scin", 9.0 * g / cm3, 5);
    material_Scin->AddElement(ele_Si, 21.9 * perCent);
    material_Scin->AddElement(ele_B, 8.8 * perCent);
    material_Scin->AddElement(ele_Na, 10.4 * perCent);
    material_Scin->AddElement(ele_Mg, 6.5 * perCent);
    material_Scin->AddElement(ele_O, 52.4 * perCent);

    if (Verbosity()) cout << "Set heavier Sciglass..." << endl;
  }
  else if ((setting > 6.) && (setting < 7.))
  {
    material_Scin = new G4Material("material_Scin", 4.5 * g / cm3, 3);
    material_Scin->AddElement(ele_Si, 21.9 * perCent);
    material_Scin->AddElement(ele_O, 52.4 * perCent);
    material_Scin->AddElement(ele_Pb, 25.7 * perCent);

    if (Verbosity()) cout << "Set Sciglass contained lead..." << endl;
  }
  else if ((setting > 7.) && (setting < 8.))
  {
    material_Scin = new G4Material("material_Scin", 4.22 * g / cm3, 4);
    material_Scin->AddElement(ele_O, 0.261);
    material_Scin->AddElement(ele_Ba, 0.3875);
    material_Scin->AddElement(ele_Si, 0.1369);
    material_Scin->AddElement(ele_Gd, 0.2146);

    if (Verbosity()) cout << "Set Sciglass from Nathaly" << endl;
  }
  else if ((setting > 8.) && (setting < 9.))
  {
    material_Scin = new G4Material("material_Scin", 3.8 * g / cm3, 3);
    material_Scin->AddElement(ele_O, 0.293);
    material_Scin->AddElement(ele_Ba, 0.502);
    material_Scin->AddElement(ele_Si, 0.205);

    if (Verbosity()) cout << "Set Sciglass from g4e" << endl;
  }

  return material_Scin;
}

//_____________________________________________________________________________
void PHG4HybridHomogeneousCalorimeterDetector::CrystalTable(G4Material* mat)
{
  //scintillation and optical properties

  //already done
  if (mat->GetMaterialPropertiesTable()) return;

  //G4cout << "OpTable::CrystalTable, " << mat->GetMaterialPropertiesTable() << G4endl;

  const G4int ntab = 2;
  G4double scin_en[] = {2.9 * eV, 3. * eV};  // 420 nm (the range is 414 - 428 nm)
  G4double scin_fast[] = {1., 1.};

  G4MaterialPropertiesTable* tab = new G4MaterialPropertiesTable();

  tab->AddProperty("FASTCOMPONENT", scin_en, scin_fast, ntab);
  tab->AddConstProperty("FASTTIMECONSTANT", 6 * ns);
  tab->AddConstProperty("SCINTILLATIONYIELD", 200 / MeV);  // 200/MEV nominal  10
  tab->AddConstProperty("RESOLUTIONSCALE", 1.);

  G4double opt_en[] = {1.551 * eV, 3.545 * eV};  // 350 - 800 nm
  G4double opt_r[] = {2.4, 2.4};
  G4double opt_abs[] = {200 * cm, 200 * cm};

  tab->AddProperty("RINDEX", opt_en, opt_r, ntab);
  tab->AddProperty("ABSLENGTH", opt_en, opt_abs, ntab);

  mat->SetMaterialPropertiesTable(tab);

}  //CrystalTable

//_____________________________________________________________________________
void PHG4HybridHomogeneousCalorimeterDetector::SurfaceTable(G4LogicalVolume* vol)
{
  //crystal optical surface

  G4OpticalSurface* surface = new G4OpticalSurface("CrystalSurface", unified, polished, dielectric_metal);
  //G4LogicalSkinSurface *csurf =
  new G4LogicalSkinSurface("CrystalSurfaceL", vol, surface);

  //surface material
  const G4int ntab = 2;
  G4double opt_en[] = {1.551 * eV, 3.545 * eV};  // 350 - 800 nm
  G4double reflectivity[] = {0.8, 0.8};
  G4double efficiency[] = {0.9, 0.9};
  G4MaterialPropertiesTable* surfmat = new G4MaterialPropertiesTable();
  surfmat->AddProperty("REFLECTIVITY", opt_en, reflectivity, ntab);
  surfmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  surface->SetMaterialPropertiesTable(surfmat);
  //csurf->DumpInfo();

}  //SurfaceTable

//_____________________________________________________________________________
void PHG4HybridHomogeneousCalorimeterDetector::MakeBoundary(G4VPhysicalVolume* crystal, G4VPhysicalVolume* opdet)
{
  //optical boundary between the crystal and optical photons detector

  G4OpticalSurface* surf = new G4OpticalSurface("OpDetS");
  //surf->SetType(dielectric_dielectric); // photons go to the detector, must have rindex defined
  surf->SetType(dielectric_metal);  // photon is absorbed when reaching the detector, no material rindex required
  //surf->SetFinish(ground);
  surf->SetFinish(polished);
  //surf->SetModel(unified);
  surf->SetModel(glisur);

  new G4LogicalBorderSurface("OpDetB", crystal, opdet, surf);

  const G4int ntab = 2;
  G4double opt_en[] = {1.551 * eV, 3.545 * eV};  // 350 - 800 nm
  //G4double reflectivity[] = {0., 0.};
  G4double reflectivity[] = {0.1, 0.1};
  //G4double reflectivity[] = {0.9, 0.9};
  //G4double reflectivity[] = {1., 1.};
  G4double efficiency[] = {1., 1.};

  G4MaterialPropertiesTable* surfmat = new G4MaterialPropertiesTable();
  surfmat->AddProperty("REFLECTIVITY", opt_en, reflectivity, ntab);
  surfmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  surf->SetMaterialPropertiesTable(surfmat);

}  //MakeBoundary

G4Material* PHG4HybridHomogeneousCalorimeterDetector::GetTedlarMaterial()
{
  static string matname = "HybridHomogeneousTedlar";
  G4Material* tedlar = G4Material::GetMaterial(matname, false);  // false suppresses warning that material does not exist
  if (!tedlar)
  {
    G4double density_tedlar = 1.43 * g / cm3;
    tedlar = new G4Material(matname, density_tedlar, 3);
    tedlar->AddElement(G4Element::GetElement("C"), 2);
    tedlar->AddElement(G4Element::GetElement("F"), 2);
    tedlar->AddElement(G4Element::GetElement("H"), 2);
  }
  return tedlar;
}

G4Material* PHG4HybridHomogeneousCalorimeterDetector::GetVM2000Material()
{
  static string matname = "HybridHomogeneousVM2000";
  G4Material* VM2000 = G4Material::GetMaterial(matname, false);  // false suppresses warning that material does not exist
  if (!VM2000)
  {
    G4double density_VM2000 = 1.43 * g / cm3;
    VM2000 = new G4Material(matname, density_VM2000, 3);
    VM2000->AddElement(G4Element::GetElement("C"), 2);
    VM2000->AddElement(G4Element::GetElement("F"), 2);
    VM2000->AddElement(G4Element::GetElement("H"), 2);

    G4MaterialPropertiesTable* mptVM2000 = new G4MaterialPropertiesTable();
    const G4int nEntriesVM2000 = 31;

    G4double photonE_VM2000[nEntriesVM2000] =
        {1.37760 * eV, 1.45864 * eV, 1.54980 * eV, 1.65312 * eV, 1.71013 * eV, 1.77120 * eV, 1.83680 * eV, 1.90745 * eV, 1.98375 * eV, 2.06640 * eV, 2.10143 * eV, 2.13766 * eV, 2.17516 * eV, 2.21400 * eV, 2.25426 * eV, 2.29600 * eV, 2.33932 * eV, 2.38431 * eV, 2.43106 * eV, 2.47968 * eV, 2.53029 * eV, 2.58300 * eV, 2.63796 * eV, 2.69531 * eV, 2.75520 * eV, 2.81782 * eV, 2.88335 * eV, 2.95200 * eV, 3.09960 * eV, 3.54241 * eV, 4.13281 * eV};
    G4double refractiveIndex_VM2000[nEntriesVM2000] =
        {1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42};
    mptVM2000->AddProperty("RINDEX", photonE_VM2000, refractiveIndex_VM2000, nEntriesVM2000);

    VM2000->SetMaterialPropertiesTable(mptVM2000);
  }
  return VM2000;
}
// <material name="PolyvinylChloride">
//   <D value="1.3" unit="g/cm3"/>
//   <fraction n="0.04838" ref="H"/>
//   <fraction n="0.38436" ref="C"/>
//   <fraction n="0.56726" ref="Cl"/>

int PHG4HybridHomogeneousCalorimeterDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  ifstream istream_mapping(m_Params->get_string_param("mappingtower"));
  if (!istream_mapping.is_open())
  {
    cout << "ERROR in PHG4HybridHomogeneousCalorimeterDetector: Failed to open mapping file " << m_Params->get_string_param("mappingtower") << endl;
    gSystem->Exit(1);
  }

  /* loop over lines in file */
  string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {
    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != string::npos)
    {
      if (Verbosity() > 0)
      {
        cout << "PHG4HybridHomogeneousCalorimeterDetector: SKIPPING line in mapping file: " << line_mapping << endl;
      }
      continue;
    }

    istringstream iss(line_mapping);

    /* If line starts with keyword Tower, add to tower positions */
    if (line_mapping.find("Tower ") != string::npos)
    {
      unsigned idx_j, idx_k, idx_l;
      G4double pos_x, pos_y, pos_z;
      G4double size_x, size_y, size_z;
      G4double rot_x, rot_y, rot_z;
      G4double dummy;
      string dummys;

      /* read string- break if error */
      if (!(iss >> dummys >> dummy >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z))
      {
        cout << "ERROR in PHG4HybridHomogeneousCalorimeterDetector: Failed to read line in mapping file " << m_Params->get_string_param("mappingtower") << endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      ostringstream towername;
      towername.str("");
      towername << _towerlogicnameprefix << "_j_" << idx_j << "_k_" << idx_k;

      /* Add Geant4 units */
      pos_x = pos_x * cm;
      pos_y = pos_y * cm;
      pos_z = pos_z * cm;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.x = pos_x;
      tower_new.y = pos_y;
      tower_new.z = pos_z;
      tower_new.idx_j = idx_j;
      tower_new.idx_k = idx_k;
      _map_tower.insert(make_pair(towername.str(), tower_new));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      string parname;
      G4double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        cout << "ERROR in PHG4HybridHomogeneousCalorimeterDetector: Failed to read line in mapping file " << m_Params->get_string_param("mappingtower") << endl;
        gSystem->Exit(1);
      }

      _map_global_parameter.insert(make_pair(parname, parval));
    }
  }

  /* Update member variables for global parameters based on parsed parameter file */
  std::map<string, G4double>::iterator parit;

  parit = _map_global_parameter.find("Gcrystal_dx");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("crystal_dx", parit->second);  // in cm
  }

  parit = _map_global_parameter.find("Gcrystal_dy");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("crystal_dy", parit->second);  // in cm
  }

  parit = _map_global_parameter.find("Gcrystal_dz");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("crystal_dz", parit->second);  // in cm
  }

  parit = _map_global_parameter.find("Gcarbon_gap");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("carbon_gap", parit->second);  // in cm
  }

  parit = _map_global_parameter.find("Gair_gap");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("air_gap", parit->second);  // in cm
  }

  parit = _map_global_parameter.find("Gr1_inner");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("rMin1", parit->second);
  }

  parit = _map_global_parameter.find("Gr1_outer");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("rMax1", parit->second);
  }

  parit = _map_global_parameter.find("Gr2_inner");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("rMin2", parit->second);
  }

  parit = _map_global_parameter.find("Gr2_outer");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("rMax2", parit->second);
  }

  parit = _map_global_parameter.find("Gdz");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("dz", parit->second);
  }

  parit = _map_global_parameter.find("Gx0");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("place_x", parit->second);
  }

  parit = _map_global_parameter.find("Gy0");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("place_y", parit->second);
  }

  parit = _map_global_parameter.find("Gz0");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("place_z", parit->second);
  }

  parit = _map_global_parameter.find("Grot_x");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("rot_x", parit->second * rad / deg);
  }

  parit = _map_global_parameter.find("Grot_y");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("rot_y", parit->second * rad / deg);
  }

  parit = _map_global_parameter.find("Grot_z");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("rot_z", parit->second * rad / deg);
  }

  parit = _map_global_parameter.find("Gmaterial");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("material", parit->second);
  }

  parit = _map_global_parameter.find("Gcolor_R");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("color_R", parit->second);
  }

  parit = _map_global_parameter.find("Gcolor_G");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("color_G", parit->second);
  }

  parit = _map_global_parameter.find("Gcolor_B");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("color_B", parit->second);
  }
  parit = _map_global_parameter.find("carbon_frame_style");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_int_param("carbon_frame_style", parit->second);
  }
  parit = _map_global_parameter.find("carbon_frame_depth");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("carbon_frame_depth", parit->second);
  }
  parit = _map_global_parameter.find("carbon_face_lip");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("carbon_face_lip", parit->second);
  }
  parit = _map_global_parameter.find("Gtedlar_thickness");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("tedlar_thickness", parit->second);
  }
  parit = _map_global_parameter.find("Greflective_foil_thickness");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("reflective_foil_thickness", parit->second);
  }
  parit = _map_global_parameter.find("sensor_dimension");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("sensor_dimension", parit->second);
  }
  parit = _map_global_parameter.find("sensor_count");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_int_param("sensor_count", parit->second);
  }
  parit = _map_global_parameter.find("sensor_thickness");
  if (parit != _map_global_parameter.end())
  {
    m_Params->set_double_param("sensor_thickness", parit->second);
  }

  return 0;
}

G4Material* PHG4HybridHomogeneousCalorimeterDetector::GetCarbonFiber()
{
  static string matname = "HybridHomogeneousCarbonFiber";
  G4Material* carbonfiber = G4Material::GetMaterial(matname, false);  // false suppresses warning that material does not exist
  if (!carbonfiber)
  {
    G4double density_carbon_fiber = 1.44 * g / cm3;
    carbonfiber = new G4Material(matname, density_carbon_fiber, 1);
    carbonfiber->AddElement(G4Element::GetElement("C"), 1);
  }
  return carbonfiber;
}
