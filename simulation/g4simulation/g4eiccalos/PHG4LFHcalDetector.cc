#include "PHG4LFHcalDetector.h"
#include "PHG4LFHcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalBorderSurface.hh>
#include <Geant4/G4LogicalSkinSurface.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PVReplica.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Torus.hh>
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

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
PHG4LFHcalDetector::PHG4LFHcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4LFHcalDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_TowerLogicNamePrefix("hHcalTower")
  , m_SuperDetector("NONE")
  , m_doLightProp(false)
{
  if (m_Params->get_string_param("mapping_file").empty())
  {
    cout << "ERROR in PHG4LFHcalDetector: No mapping file specified. Abort detector construction." << endl;
    cout << "Please run set_string_param(\"mapping_file\", std::string filename ) first." << endl;
    gSystem->Exit(1);
  }
  /* Read parameters for detector construction and mappign from file */
  ParseParametersFromTable();
}
//_______________________________________________________________________
int PHG4LFHcalDetector::IsInLFHcal(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* mylogvol = volume->GetLogicalVolume();
  if (m_ActiveFlag)
  {
    if (m_ScintiLogicalVolSet.find(mylogvol) != m_ScintiLogicalVolSet.end())
    {
      return 1;
    }
  }

  if (m_AbsorberActiveFlag)
  {
    if (m_AbsorberLogicalVolSet.find(mylogvol) != m_AbsorberLogicalVolSet.end())
    {
      return -1;
    }
  }
  return 0;
}

//_______________________________________________________________________
void PHG4LFHcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4LFHcalDetector: Begin Construction" << std::endl;
  }

  // /* Read parameters for detector construction and mappign from file */
  // ParseParametersFromTable();

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));

  G4VSolid* beampipe_cutout = new G4Cons("LFHCAL_beampipe_cutout",
                                         0, m_Params->get_double_param("rMin1") * cm,
                                         0, m_Params->get_double_param("rMin2") * cm,
                                         m_Params->get_double_param("dz") * cm / 2.0,
                                         0, 2 * M_PI);
  G4VSolid* hcal_envelope_solid = new G4Cons("LFHCAL_envelope_solid_cutout",
                                             0, m_Params->get_double_param("rMax1") * cm,
                                             0, m_Params->get_double_param("rMax2") * cm,
                                             m_Params->get_double_param("dz") * cm / 2.0,
                                             0, 2 * M_PI);
  hcal_envelope_solid = new G4SubtractionSolid(G4String("LFHCAL_envelope_solid"), hcal_envelope_solid, beampipe_cutout, 0, G4ThreeVector(m_Params->get_double_param("xoffset") * cm, m_Params->get_double_param("yoffset") * cm, 0.));

  G4LogicalVolume* hcal_envelope_log = new G4LogicalVolume(hcal_envelope_solid, WorldMaterial, "hLFHCAL_envelope", 0, 0, 0);

  m_DisplayAction->AddVolume(hcal_envelope_log, "LFHcalEnvelope");

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  hcal_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  hcal_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);

  /* Place envelope cone in simulation */
  string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm,
                                                           m_Params->get_double_param("place_y") * cm,
                                                           m_Params->get_double_param("place_z") * cm)),
                    hcal_envelope_log, name_envelope, logicWorld, 0, false, OverlapCheck());

  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  /* Place calorimeter tower within envelope */
  PlaceTower(hcal_envelope_log, singletower);

  return;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4LFHcalDetector::ConstructTower()
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4LFHcalDetector: Build logical volume for single tower..." << std::endl;
  }

  //**********************************************************************************************
  /* read variables from external inputs */
  //**********************************************************************************************
  G4double TowerDx = m_Params->get_double_param("tower_dx") * cm;
  G4double TowerDy = m_Params->get_double_param("tower_dy") * cm;
  G4double TowerDz = m_Params->get_double_param("tower_dz") * cm;
  G4double WlsDw = m_Params->get_double_param("wls_dw") * cm;
  G4double thick_frame_width = m_Params->get_double_param("frame_width") * cm;
  G4double thin_frame_width = 0;
  if (thick_frame_width > 0)
    thin_frame_width = 0.1 * mm;
  // double width_coating    = m_Params->get_double_param("width_coating") * cm;
  G4double thickness_absorber = m_Params->get_double_param("thickness_absorber") * cm;
  G4double thickness_scintillator = m_Params->get_double_param("thickness_scintillator") * cm;
  G4int nlayers = TowerDz / (thickness_absorber + thickness_scintillator);
  G4Material* material_scintillator = GetScintillatorMaterial();  //G4Material::GetMaterial(m_Params->get_string_param("scintillator"));
  G4Material* material_absorber = GetDetectorMaterial(m_Params->get_string_param("absorber"));
  G4Material* material_absorber_W = GetDetectorMaterial(m_Params->get_string_param("absorber_W"));
  G4Material* material_wls = GetWLSFiberMaterial();  //GetDetectorMaterial(m_Params->get_string_param("scintillator"));
  int embed_fiber = m_Params->get_int_param("embed_fiber");
  G4double fiber_thickness = 0.2 * mm;
  // G4double fiber_thickness        = 1.0*mm;
  //**********************************************************************************************
  /* create logical volume for single tower */
  //**********************************************************************************************
  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));
  G4VSolid* single_tower_solid = new G4Box("single_tower_solid",
                                           TowerDx / 2.0,
                                           TowerDy / 2.0,
                                           TowerDz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            WorldMaterial,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  G4double notch_length = 2.0 * cm;
  std::vector<G4TwoVector> poligon = {
      {(-TowerDx) / 2.0 + thin_frame_width, (-TowerDy) / 2.0 + thick_frame_width},
      {(-TowerDx) / 2.0 + thin_frame_width, (TowerDy) / 2.0 - thick_frame_width},
      {(TowerDx - WlsDw) / 2.0 - thin_frame_width, (TowerDy) / 2.0 - thick_frame_width},
      {(TowerDx - WlsDw) / 2.0 - thin_frame_width, (-TowerDy + notch_length) / 2.0 + thick_frame_width},
      {(TowerDx) / 2.0 - thin_frame_width, (-TowerDy + notch_length) / 2.0 + thick_frame_width},
      {(TowerDx) / 2.0 - thin_frame_width, (-TowerDy) / 2.0 + thick_frame_width}};
  //**********************************************************************************************
  /* create logical and geometry volumes for minitower read-out unit */
  //**********************************************************************************************
  std::vector<G4ExtrudedSolid::ZSection> zsections_miniblock = {{-(thickness_absorber + thickness_scintillator) / 2.0, {0, 0}, 1.0}, {(thickness_absorber + thickness_scintillator) / 2.0, {0, 0}, 1.}};
  G4VSolid* miniblock_solid = new G4ExtrudedSolid("miniblock_solid", poligon, zsections_miniblock);

  G4LogicalVolume* miniblock_logic = new G4LogicalVolume(miniblock_solid,
                                                         WorldMaterial,
                                                         "miniblock_logic",
                                                         0, 0, 0);
  G4LogicalVolume* miniblock_W_logic = new G4LogicalVolume(miniblock_solid,
                                                           WorldMaterial,
                                                           "miniblock_W_logic",
                                                           0, 0, 0);
  m_DisplayAction->AddVolume(miniblock_logic, "Invisible");
  m_DisplayAction->AddVolume(miniblock_W_logic, "Invisible");
  //**********************************************************************************************
  /* create logical & geometry volumes for scintillator and absorber plates to place inside mini read-out unit */
  //**********************************************************************************************
  std::vector<G4ExtrudedSolid::ZSection> zsections_scintillator = {{-thickness_scintillator / 2, {0, 0}, 1.0}, {thickness_scintillator / 2, {0, 0}, 1.}};

  G4VSolid* solid_scintillator = new G4ExtrudedSolid("single_plate_scintillator", poligon, zsections_scintillator);

  G4VPhysicalVolume* physvol_fiber_loop_0 = nullptr;
  G4VPhysicalVolume* physvol_fiber_loop_1 = nullptr;
  G4VPhysicalVolume* physvol_fiber_straight_1 = nullptr;
  G4LogicalVolume* logic_embed_fiber_loop = nullptr;
  G4LogicalVolume* logic_embed_fiber_loop_2 = nullptr;
  G4LogicalVolume* logic_embed_fiber_straight_1 = nullptr;
  G4double cutout_margin = 0.1 * fiber_thickness;
  G4double fiber_bending_R_1 = notch_length / 2.0 - (TowerDy / 2 - thick_frame_width - (0.9 * (TowerDx - WlsDw - 2 * thin_frame_width) / 2.0) - (2 * cutout_margin));
  G4double fiber_bending_R = 1.0 * cm;
  if (embed_fiber)
  {
    G4VSolid* solid_embed_fiber_loop = new G4Torus("solid_embed_fiber_loop",
                                                   0, fiber_thickness / 2.0,
                                                   (0.9 * (TowerDx - WlsDw - 2 * thin_frame_width) / 2.0) - (fiber_thickness + 2 * cutout_margin) / 2.0,
                                                   -0.4 * M_PI, 1.9 * M_PI);
    G4VSolid* solid_embed_fiber_loop_2 = new G4Torus("solid_embed_fiber_loop_2",
                                                     0, fiber_thickness / 2.0,
                                                     (2 * fiber_bending_R_1) / 2.0 - (fiber_thickness + 2 * cutout_margin) / 2.0,
                                                     -0.5 * M_PI, 0.5 * M_PI);
    G4VSolid* solid_embed_fiber_straight_1 = new G4Tubs("solid_embed_fiber_straight_1",
                                                        0,
                                                        fiber_thickness / 2.0,
                                                        // (TowerDx/2.0 - fiber_bending_R_1 - thin_frame_width) / 2.0,
                                                        ((TowerDx - WlsDw + fiber_thickness) / 2.0 - (fiber_bending_R_1)) / 2.0,
                                                        0, 2.0 * M_PI);
    solid_scintillator = new G4SubtractionSolid(G4String("single_plate_scintillator_cu"),
                                                solid_scintillator, solid_embed_fiber_loop,
                                                0, G4ThreeVector(0, 0, -thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0));
    solid_scintillator = new G4SubtractionSolid(G4String("single_plate_scintillator_cu2"),
                                                solid_scintillator, solid_embed_fiber_loop_2,
                                                0, G4ThreeVector((TowerDx - WlsDw + fiber_thickness / 2.0) / 2.0 - (fiber_bending_R_1) + (thin_frame_width) + cutout_margin, -0.9 * (TowerDx - WlsDw - thin_frame_width) / 2.0 + (fiber_bending_R_1), -thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0));
    G4RotationMatrix* wls_rotm_fibrcu = new G4RotationMatrix();
    wls_rotm_fibrcu->rotateY(90 * deg);
    solid_scintillator = new G4SubtractionSolid(G4String("single_plate_scintillator_cu3"),
                                                solid_scintillator, solid_embed_fiber_straight_1,
                                                wls_rotm_fibrcu, G4ThreeVector(((TowerDx - WlsDw + fiber_thickness) / 2.0 - (fiber_bending_R_1)) / 2.0,  //(TowerDx/2-((TowerDx - WlsDw+fiber_thickness/2.0)/2.0-(fiber_bending_R_1)+(thin_frame_width)+cutout_margin))/2+fiber_bending_R_1/2+(fiber_thickness),
                                                                               -0.9 * (TowerDx - WlsDw - thin_frame_width) / 2.0 + (fiber_thickness + 2 * cutout_margin) / 2.0, -thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0));
    logic_embed_fiber_loop = new G4LogicalVolume(solid_embed_fiber_loop,
                                                 material_wls, "logic_embed_fiber_loop",
                                                 0, 0, 0);
    logic_embed_fiber_loop_2 = new G4LogicalVolume(solid_embed_fiber_loop_2,
                                                   material_wls, "logic_embed_fiber_loop_2",
                                                   0, 0, 0);
    logic_embed_fiber_straight_1 = new G4LogicalVolume(solid_embed_fiber_straight_1,
                                                       material_wls, "logic_embed_fiber_straight_1",
                                                       0, 0, 0);
    // if(m_doLightProp){
    //   SurfaceTable(logic_embed_fiber_loop);
    //   SurfaceTable(logic_embed_fiber_loop_2);
    //   SurfaceTable(logic_embed_fiber_straight_1);
    // }
    physvol_fiber_loop_0 = new G4PVPlacement(0, G4ThreeVector(0, 0, (thickness_absorber) / 2. - thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0),
                                             logic_embed_fiber_loop,
                                             "embed_fiber_loop_placed",
                                             miniblock_logic,
                                             0, 0, OverlapCheck());
    physvol_fiber_loop_1 = new G4PVPlacement(0, G4ThreeVector(
                                                    // (TowerDx - WlsDw+fiber_thickness/2.0)/2.0-(fiber_bending_R_1),
                                                    (TowerDx - WlsDw + fiber_thickness / 2.0) / 2.0 - (fiber_bending_R_1) + (thin_frame_width) + cutout_margin,
                                                    // (TowerDx)/2.0-(fiber_bending_R_1) - thin_frame_width,
                                                    -0.9 * (TowerDx - WlsDw - thin_frame_width) / 2.0 + (fiber_bending_R_1), (thickness_absorber) / 2. - thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0),
                                             logic_embed_fiber_loop_2,
                                             "embed_fiber_loop_2_placed",
                                             miniblock_logic,
                                             0, 0, OverlapCheck());
    G4RotationMatrix* wls_rotm_fibr = new G4RotationMatrix();
    wls_rotm_fibr->rotateY(90 * deg);
    physvol_fiber_straight_1 = new G4PVPlacement(wls_rotm_fibr, G4ThreeVector(
                                                                    // (TowerDx/2.0 - fiber_bending_R_1 - thin_frame_width) / 2.0,
                                                                    ((TowerDx - WlsDw + fiber_thickness) / 2.0 - (fiber_bending_R_1)) / 2.0,
                                                                    // (TowerDx/2-((TowerDx - WlsDw+fiber_thickness/2.0)/2.0-(fiber_bending_R_1)+(thin_frame_width)+cutout_margin))/2,
                                                                    -0.9 * (TowerDx - WlsDw - thin_frame_width) / 2.0 + (fiber_thickness + 2 * cutout_margin) / 2.0, (thickness_absorber) / 2. - thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0),
                                                 logic_embed_fiber_straight_1,
                                                 "embed_fiber_straight_1_placed",
                                                 miniblock_logic,
                                                 0, 0, OverlapCheck());
    if (m_doLightProp)
    {
      MakeBoundary(physvol_fiber_loop_0, physvol_fiber_loop_1, true);
      MakeBoundary(physvol_fiber_loop_1, physvol_fiber_straight_1, true);
    }
    m_DisplayAction->AddVolume(logic_embed_fiber_loop, "WLSfiber");
    m_DisplayAction->AddVolume(logic_embed_fiber_loop_2, "WLSfiber");
    m_DisplayAction->AddVolume(logic_embed_fiber_straight_1, "WLSfiber");
  }
  // G4VSolid* solid_absorber = new G4Box("single_plate_absorber_solid",
  //                                      (TowerDx - WlsDw) / 2.0,
  //                                      (TowerDy ) / 2.0,
  //                                      thickness_absorber / 2.0);

  std::vector<G4ExtrudedSolid::ZSection> zsections_absorber = {{-thickness_absorber / 2, {0, 0}, 1.0}, {thickness_absorber / 2, {0, 0}, 1.}};
  G4VSolid* solid_absorber = new G4ExtrudedSolid("single_plate_absorber_solid", poligon, zsections_absorber);
  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "single_plate_absorber_logic",
                                                        0, 0, 0);
  m_AbsorberLogicalVolSet.insert(logic_absorber);
  G4LogicalVolume* logic_absorber_W = new G4LogicalVolume(solid_absorber,
                                                          material_absorber_W,
                                                          "single_plate_absorber_W_logic",
                                                          0, 0, 0);
  m_AbsorberLogicalVolSet.insert(logic_absorber_W);
  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                     material_scintillator,
                                                     "hLFHCAL_scintillator_plate_logic",
                                                     0, 0, 0);
  m_ScintiLogicalVolSet.insert(logic_scint);
  // if(m_doLightProp){
  //   SurfaceTable(logic_scint);
  // }
  m_DisplayAction->AddVolume(logic_absorber, "Absorber");
  m_DisplayAction->AddVolume(logic_absorber_W, "Absorber_W");
  m_DisplayAction->AddVolume(logic_scint, "Scintillator");
  string name_absorber = m_TowerLogicNamePrefix + "_single_plate_absorber";
  string name_absorber_W = m_TowerLogicNamePrefix + "_single_plate_absorber_W";
  string name_scintillator = m_TowerLogicNamePrefix + "_single_plate_scintillator";
  // place Steel absorber and scintillator in miniblock
  new G4PVPlacement(0, G4ThreeVector(0, 0, -thickness_scintillator / 2),
                    logic_absorber,
                    name_absorber,
                    miniblock_logic,
                    0, 0, OverlapCheck());

  G4VPhysicalVolume* scintillator_placed = new G4PVPlacement(0, G4ThreeVector(0, 0, (thickness_absorber) / 2.),
                                                             logic_scint,
                                                             name_scintillator,
                                                             miniblock_logic,
                                                             0, 0, OverlapCheck());
  if (m_doLightProp && embed_fiber)
  {
    MakeBoundary_Fiber_Scint(physvol_fiber_loop_0, scintillator_placed);
    MakeBoundary_Fiber_Scint(physvol_fiber_loop_1, scintillator_placed);
    MakeBoundary_Fiber_Scint(physvol_fiber_straight_1, scintillator_placed);
  }
  // place Tungsten absorber and scintillator in a separate miniblock
  new G4PVPlacement(0, G4ThreeVector(0, 0, -thickness_scintillator / 2),
                    logic_absorber_W,
                    name_absorber,
                    miniblock_W_logic,
                    0, 0, OverlapCheck());

  new G4PVPlacement(0, G4ThreeVector(0, 0, (thickness_absorber) / 2.),
                    logic_scint,
                    name_scintillator,
                    miniblock_W_logic,
                    0, 0, OverlapCheck());
  if (embed_fiber)
  {
    G4VPhysicalVolume* physvol_fiber_loop_W_0 = new G4PVPlacement(0, G4ThreeVector(0, 0, (thickness_absorber) / 2. - thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0),
                                                                  logic_embed_fiber_loop,
                                                                  "embed_fiber_loop_W_placed",
                                                                  miniblock_W_logic,
                                                                  0, 0, OverlapCheck());
    G4VPhysicalVolume* physvol_fiber_loop_W_1 = new G4PVPlacement(0, G4ThreeVector((TowerDx - WlsDw + fiber_thickness / 2.0) / 2.0 - (fiber_bending_R_1) + (thin_frame_width) + cutout_margin,
                                                                                   // (TowerDx)/2.0-(fiber_bending_R_1) - thin_frame_width,
                                                                                   -0.9 * (TowerDx - WlsDw - thin_frame_width) / 2.0 + (fiber_bending_R_1), (thickness_absorber) / 2. - thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0),
                                                                  logic_embed_fiber_loop_2,
                                                                  "embed_fiber_loop_W_2_placed",
                                                                  miniblock_W_logic,
                                                                  0, 0, OverlapCheck());
    G4RotationMatrix* wls_rotm_fibr = new G4RotationMatrix();
    wls_rotm_fibr->rotateY(90 * deg);
    G4VPhysicalVolume* physvol_fiber_loop_W_2 = new G4PVPlacement(wls_rotm_fibr, G4ThreeVector(((TowerDx - WlsDw + fiber_thickness) / 2.0 - (fiber_bending_R_1)) / 2.0,
                                                                                               // (TowerDx/2-((TowerDx - WlsDw+fiber_thickness/2.0)/2.0-(fiber_bending_R_1)+(thin_frame_width)+cutout_margin))/2,
                                                                                               -0.9 * (TowerDx - WlsDw - thin_frame_width) / 2.0 + (fiber_thickness + 2 * cutout_margin) / 2.0, (thickness_absorber) / 2. - thickness_scintillator / 2 + (fiber_thickness + cutout_margin) / 2.0),
                                                                  logic_embed_fiber_straight_1,
                                                                  "embed_fiber_straight_W_1_placed",
                                                                  miniblock_W_logic,
                                                                  0, 0, OverlapCheck());
    if (m_doLightProp)
    {
      MakeBoundary(physvol_fiber_loop_W_0, physvol_fiber_loop_W_1, true);
      MakeBoundary(physvol_fiber_loop_W_1, physvol_fiber_loop_W_2, true);
    }
  }
  //**********************************************************************************************
  /* create wavelength shifting "fiber-block"  */
  //**********************************************************************************************
  if (thin_frame_width > 0)
  {
    G4Material* material_frame = GetDetectorMaterial("G4_Fe");

    G4VSolid* solid_frame_plate = new G4Box("single_plate_frame",
                                            (thin_frame_width) / 2, (TowerDy - 2 * thick_frame_width) / 2, TowerDz / 2);
    G4LogicalVolume* logic_frame = new G4LogicalVolume(solid_frame_plate,
                                                       material_frame,
                                                       "hLFHCAL_frame_plate_logic",
                                                       0, 0, 0);
    m_DisplayAction->AddVolume(logic_frame, "Frame");

    new G4PVPlacement(0, G4ThreeVector(-TowerDx / 2 + thin_frame_width / 2, 0, 0.),
                      logic_frame,
                      m_TowerLogicNamePrefix + "_single_plate_frame_left",
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector(TowerDx / 2 - thin_frame_width / 2, 0, 0.),
                      logic_frame,
                      m_TowerLogicNamePrefix + "_single_plate_frame_right",
                      single_tower_logic,
                      0, 0, OverlapCheck());

    G4VSolid* solid_cover_plate = new G4Box("single_plate_cover",
                                            TowerDx / 2, thick_frame_width / 2, TowerDz / 2);
    G4LogicalVolume* logic_cover = new G4LogicalVolume(solid_cover_plate,
                                                       material_frame,
                                                       "hLFHCAL_cover_plate_logic",
                                                       0, 0, 0);
    m_DisplayAction->AddVolume(logic_cover, "Frame");

    new G4PVPlacement(0, G4ThreeVector(0, -TowerDy / 2 + thick_frame_width / 2, 0.),
                      logic_cover,
                      m_TowerLogicNamePrefix + "_single_plate_cover_top",
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector(0, TowerDy / 2 - thick_frame_width / 2, 0.),
                      logic_cover,
                      m_TowerLogicNamePrefix + "_single_plate_cover_bottom",
                      single_tower_logic,
                      0, 0, OverlapCheck());
  }

  //**********************************************************************************************
  /* create logical volume for single tower */
  //**********************************************************************************************
  double SteelTowerLength = TowerDz;
  double WTowerLength = 0.;
  int nLayersSteel = nlayers;
  int nLayersTungsten = 0;
  if (m_Params->get_int_param("usetailcatcher"))
  {
    SteelTowerLength -= 10 * (thickness_absorber + thickness_scintillator);
    nLayersSteel -= 10;
    nLayersTungsten = 10;
    WTowerLength = 10 * (thickness_absorber + thickness_scintillator);
    std::cout << "using 10 layer tungsten tailcatcher in LFHCAL" << std::endl;
  }

  std::vector<G4ExtrudedSolid::ZSection> zsections_steeltower = {{-(SteelTowerLength) / 2.0, {0, 0}, 1.0}, {(SteelTowerLength) / 2.0, {0, 0}, 1.}};
  G4VSolid* single_tower_solidRep = new G4ExtrudedSolid("single_tower_solidRep", poligon, zsections_steeltower);

  G4LogicalVolume* single_tower_logicRep = new G4LogicalVolume(single_tower_solidRep,
                                                               WorldMaterial,
                                                               "single_tower_logicRep",
                                                               0, 0, 0);
  string name_tower = m_TowerLogicNamePrefix;
  new G4PVReplica(name_tower, miniblock_logic, single_tower_logicRep,
                  kZAxis, nLayersSteel, thickness_absorber + thickness_scintillator, 0);

  new G4PVPlacement(0, G4ThreeVector(0, 0, -WTowerLength / 2),
                    single_tower_logicRep,
                    name_tower,
                    single_tower_logic,
                    0, 0, OverlapCheck());

  m_DisplayAction->AddVolume(single_tower_logicRep, "SingleTower");
  m_DisplayAction->AddVolume(single_tower_logic, "SingleTower");

  if (m_Params->get_int_param("usetailcatcher"))
  {
    std::vector<G4ExtrudedSolid::ZSection> zsections_wtower = {{-(WTowerLength) / 2.0, {0, 0}, 1.0}, {(WTowerLength) / 2.0, {0, 0}, 1.}};
    G4VSolid* single_tower_W_solidRep = new G4ExtrudedSolid("single_tower_W_solidRep", poligon, zsections_wtower);

    G4LogicalVolume* single_tower_W_logicRep = new G4LogicalVolume(single_tower_W_solidRep,
                                                                   WorldMaterial,
                                                                   "single_tower_W_logicRep",
                                                                   0, 0, 0);
    string name_tower_W = m_TowerLogicNamePrefix + "_W";
    new G4PVReplica(name_tower_W, miniblock_W_logic, single_tower_W_logicRep,
                    kZAxis, nLayersTungsten, thickness_absorber + thickness_scintillator, 0);

    new G4PVPlacement(0, G4ThreeVector(0, 0, SteelTowerLength / 2),
                      single_tower_W_logicRep,
                      name_tower_W,
                      single_tower_logic,
                      0, 0, OverlapCheck());

    m_DisplayAction->AddVolume(single_tower_W_logicRep, "SingleTower_W");
  }
  if (embed_fiber)
  {
    G4double spacer_width = 1.5 * mm;
    G4double add_spacing = spacer_width + 0.1 * mm;
    G4double extraspacing = add_spacing;
    for (int ilay = 0; ilay < nlayers; ilay++)
    {
      if ((SteelTowerLength + WTowerLength - ilay * (thickness_absorber + thickness_scintillator) - thickness_absorber - thickness_scintillator / 2 - fiber_bending_R / 2) / 2.0 > 0)
      {
        G4VSolid* solid_long_fiber_tmp = new G4Tubs("solid_long_fiber_tmp_" + std::to_string(ilay),
                                                    0,
                                                    fiber_thickness / 2.0,
                                                    (SteelTowerLength + WTowerLength - ilay * (thickness_absorber + thickness_scintillator) - thickness_absorber - fiber_thickness / 2 - fiber_bending_R / 2) / 2.0,
                                                    0, 2.0 * M_PI);

        G4LogicalVolume* logic_long_fiber_tmp = new G4LogicalVolume(solid_long_fiber_tmp,
                                                                    material_wls,
                                                                    "logic_long_fiber_tmp" + std::to_string(ilay),
                                                                    0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector((TowerDx - WlsDw + fiber_thickness / 2.0) / 2.0, TowerDy / 2 - thick_frame_width - fiber_thickness - ilay * (1.01 * fiber_thickness) - extraspacing, (ilay * (thickness_absorber + thickness_scintillator) + thickness_absorber + fiber_thickness / 2 + fiber_bending_R / 2) / 2.0),
                          logic_long_fiber_tmp,
                          "phys_long_fiber_tmp" + std::to_string(ilay),
                          single_tower_logic,
                          0, 0, OverlapCheck());
        m_DisplayAction->AddVolume(logic_long_fiber_tmp, "WLSfiber");
      }

      // short fiber piece from loop to backward readout
      G4double shortup_fiber_length = TowerDy / 2 - thick_frame_width - fiber_thickness - ilay * (1.01 * fiber_thickness) - extraspacing - fiber_bending_R / 2 + TowerDy / 2 - thick_frame_width - notch_length / 2;
      G4VSolid* solid_shortup_fiber_tmp = new G4Tubs("solid_shortup_fiber_tmp_" + std::to_string(ilay),
                                                     0,
                                                     fiber_thickness / 2.0,
                                                     shortup_fiber_length / 2.0,
                                                     0, 2.0 * M_PI);
      G4LogicalVolume* logic_shortup_fiber_tmp = new G4LogicalVolume(solid_shortup_fiber_tmp,
                                                                     material_wls,
                                                                     "logic_shortup_fiber_tmp" + std::to_string(ilay),
                                                                     0, 0, 0);

      G4RotationMatrix* wls_rotm_fibr3 = new G4RotationMatrix();
      wls_rotm_fibr3->rotateX(90 * deg);
      new G4PVPlacement(wls_rotm_fibr3, G4ThreeVector((TowerDx - WlsDw + fiber_thickness / 2.0) / 2.0, TowerDy / 2 - thick_frame_width - fiber_thickness - ilay * (1.01 * fiber_thickness) - extraspacing - fiber_bending_R / 2 - (shortup_fiber_length / 2.0), -TowerDz / 2 + ilay * (thickness_absorber + thickness_scintillator) + thickness_absorber + fiber_thickness / 2),
                        logic_shortup_fiber_tmp,
                        "phys_shortup_fiber_tmp" + std::to_string(ilay),
                        single_tower_logic,
                        0, 0, OverlapCheck());
      m_DisplayAction->AddVolume(logic_shortup_fiber_tmp, "WLSfiber");

      // curved fiber piece to backward readout
      G4VSolid* solid_fiber_curved = new G4Torus("solid_fiber_curved",
                                                 0, fiber_thickness / 2.0,
                                                 (fiber_bending_R) / 2.0,
                                                 0.5 * M_PI, 0.5 * M_PI);
      if (ilay == nlayers - 1)
      {
        solid_fiber_curved = new G4Torus("solid_fiber_curved_lastlayer",
                                         0, fiber_thickness / 2.0,
                                         (fiber_bending_R) / 2.0,
                                         0.6 * M_PI, 0.4 * M_PI);
      }
      G4LogicalVolume* logic_fiber_curved = new G4LogicalVolume(solid_fiber_curved,
                                                                material_wls,
                                                                "logic_fiber_curved" + std::to_string(ilay),
                                                                0, 0, 0);

      G4RotationMatrix* wls_rotm_fibr2 = new G4RotationMatrix();
      wls_rotm_fibr2->rotateY(90 * deg);
      new G4PVPlacement(wls_rotm_fibr2, G4ThreeVector((TowerDx - WlsDw + fiber_thickness / 2.0) / 2.0, TowerDy / 2 - thick_frame_width - fiber_thickness - ilay * (1.01 * fiber_thickness) - extraspacing - fiber_bending_R / 2, -TowerDz / 2 + ilay * (thickness_absorber + thickness_scintillator) + thickness_absorber + fiber_thickness / 2 + fiber_bending_R / 2),
                        logic_fiber_curved,
                        "phys_fiber_curved" + std::to_string(ilay),
                        single_tower_logic,
                        0, 0, OverlapCheck());
      m_DisplayAction->AddVolume(logic_fiber_curved, "WLSfiber");
      if ((ilay + 1) % 10 == 0)
      {
        if (((SteelTowerLength + WTowerLength - ilay * (thickness_absorber + thickness_scintillator) - thickness_absorber - thickness_scintillator / 2 - thickness_absorber / 2) / 2.0) > 0)
        {
          G4VSolid* solid_spacer_tmp = new G4Box("solid_spacer_tmp_" + std::to_string(ilay),
                                                 (WlsDw) / 4,
                                                 spacer_width / 2.0,
                                                 (SteelTowerLength + WTowerLength - ilay * (thickness_absorber + thickness_scintillator) - thickness_absorber - thickness_scintillator / 2 - thickness_absorber / 2) / 2.0);

          G4LogicalVolume* logic_spacer_tmp = new G4LogicalVolume(solid_spacer_tmp,
                                                                  G4Material::GetMaterial("CFRP_INTT"),  // carbon fiber + epoxy
                                                                  "logic_spacer_tmp" + std::to_string(ilay),
                                                                  0, 0, 0);
          new G4PVPlacement(0, G4ThreeVector((TowerDx - WlsDw / 2) / 2.0 - thin_frame_width, TowerDy / 2 - thick_frame_width - fiber_thickness - (ilay + 0.5) * (1.01 * fiber_thickness) - extraspacing - add_spacing / 2, (ilay * (thickness_absorber + thickness_scintillator) + thickness_absorber + thickness_scintillator / 2 + thickness_absorber / 2) / 2.0),
                            logic_spacer_tmp,
                            "phys_spacer_tmp" + std::to_string(ilay),
                            single_tower_logic,
                            0, 0, OverlapCheck());
          m_DisplayAction->AddVolume(logic_spacer_tmp, "Spacer");
        }
        extraspacing += add_spacing;
      }
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << "PHG4LFHcalDetector: Building logical volume for single tower done." << std::endl;
  }

  return single_tower_logic;
}

int PHG4LFHcalDetector::PlaceTower(G4LogicalVolume* hcalenvelope, G4LogicalVolume* singletower)
{
  /* Loop over all tower positions in vector and place tower */
  for (std::map<std::string, towerposition>::iterator iterator = m_TowerPostionMap.begin(); iterator != m_TowerPostionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4LFHcalDetector: Place tower " << iterator->first
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k
                << " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << std::endl;
    }

    int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;
    //     std::cout << "tower " << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k << " cpNo: " << copyno << std::endl;
    new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                      singletower,
                      iterator->first,
                      hcalenvelope,
                      0, copyno, OverlapCheck());
  }

  return 0;
}

//_______________________________________________________________________
G4Material* PHG4LFHcalDetector::GetScintillatorMaterial()
{
  G4double density;
  G4int ncomponents;
  G4Material* material_ScintFEMC = new G4Material("PolystyreneFEMC", density = 1.03 * g / cm3, ncomponents = 2);
  material_ScintFEMC->AddElement(G4Element::GetElement("C"), 8);
  material_ScintFEMC->AddElement(G4Element::GetElement("H"), 8);

  if (m_doLightProp)
  {
    // const G4int ntab = 31;
    // G4double opt_en[] =
    //   { 1.37760*eV, 1.45864*eV, 1.54980*eV, 1.65312*eV, 1.71013*eV, 1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV,
    //     2.10143*eV, 2.13766*eV, 2.17516*eV, 2.21400*eV, 2.25426*eV, 2.29600*eV, 2.33932*eV, 2.38431*eV, 2.43106*eV, 2.47968*eV,
    //     2.53029*eV, 2.58300*eV, 2.63796*eV, 2.69531*eV, 2.75520*eV, 2.81782*eV, 2.88335*eV, 2.95200*eV, 3.09960*eV, 3.54241*eV,
    //     4.13281*eV }; // 350 - 800 nm
    // G4double opt_abs[] =
    //   { 2.714*m, 3.619*m, 5.791*m, 4.343*m, 7.896*m, 5.429*m, 36.19*m, 17.37*m, 36.19*m, 5.429*m,
    //     13.00*m, 14.50*m, 16.00*m, 18.00*m, 16.50*m, 17.00*m, 14.00*m, 16.00*m, 15.00*m, 14.50*m,
    //     13.00*m, 12.00*m, 10.00*m, 8.000*m, 7.238*m, 4.000*m, 1.200*m, 0.500*m, 0.200*m, 0.200*m,
    //     0.100*m };

    const G4int nEntries = 50;
    G4double photonEnergy[nEntries] =
        {2.00 * eV, 2.03 * eV, 2.06 * eV, 2.09 * eV, 2.12 * eV,
         2.15 * eV, 2.18 * eV, 2.21 * eV, 2.24 * eV, 2.27 * eV,
         2.30 * eV, 2.33 * eV, 2.36 * eV, 2.39 * eV, 2.42 * eV,
         2.45 * eV, 2.48 * eV, 2.51 * eV, 2.54 * eV, 2.57 * eV,
         2.60 * eV, 2.63 * eV, 2.66 * eV, 2.69 * eV, 2.72 * eV,
         2.75 * eV, 2.78 * eV, 2.81 * eV, 2.84 * eV, 2.87 * eV,
         2.90 * eV, 2.93 * eV, 2.96 * eV, 2.99 * eV, 3.02 * eV,
         3.05 * eV, 3.08 * eV, 3.11 * eV, 3.14 * eV, 3.17 * eV,
         3.20 * eV, 3.23 * eV, 3.26 * eV, 3.29 * eV, 3.32 * eV,
         3.35 * eV, 3.38 * eV, 3.41 * eV, 3.44 * eV, 3.47 * eV};
    G4double scintilFast[nEntries] =
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    const G4int ntab = 4;

    G4double wls_Energy[] = {2.00 * eV, 2.87 * eV, 2.90 * eV,
                             3.47 * eV};

    G4double rIndexPstyrene[] = {1.5, 1.5, 1.5, 1.5};
    G4double absorption1[] = {2. * cm, 2. * cm, 2. * cm, 2. * cm};
    // G4double scintilFast[]    = { 0.0, 0.0, 1.0, 1.0 };
    G4MaterialPropertiesTable* fMPTPStyrene = new G4MaterialPropertiesTable();
    fMPTPStyrene->AddProperty("RINDEX", wls_Energy, rIndexPstyrene, ntab);
    fMPTPStyrene->AddProperty("ABSLENGTH", wls_Energy, absorption1, ntab);
    // fMPTPStyrene->AddProperty("ABSLENGTH", opt_en, opt_abs, ntab);
    // fMPTPStyrene->AddProperty("SCINTILLATIONCOMPONENT1", wls_Energy, scintilFast,ntab);
    fMPTPStyrene->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintilFast, nEntries);
    fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
    fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fMPTPStyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT", 10. * ns);
    // fMPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);

    // fMPTPStyrene->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,nEntries);
    material_ScintFEMC->SetMaterialPropertiesTable(fMPTPStyrene);
  }
  // Set the Birks Constant for the Polystyrene scintillator
  material_ScintFEMC->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  return material_ScintFEMC;
}

//_______________________________________________________________________
G4Material* PHG4LFHcalDetector::GetCoatingMaterial()
{
  //--------------------------------------------------
  // TiO2
  //--------------------------------------------------
  G4double density, fractionmass;
  G4double a, z;
  G4int ncomponents;
  G4Material* material_TiO2 = new G4Material("TiO2_FEMC", density = 1.52 * g / cm3, ncomponents = 2);
  G4Element* Ti = new G4Element("Titanium", "Ti", z=22, a=  47.8670*g/mole);
  material_TiO2->AddElement(Ti, 1);
  // material_TiO2->AddElement(G4Element::GetElement("Ti"), 1);
  material_TiO2->AddElement(G4Element::GetElement("O"), 2);

  //--------------------------------------------------
  // Scintillator Coating - 15% TiO2 and 85% polystyrene by weight.
  //--------------------------------------------------
  //Coating_FEMC (Glass + Epoxy)
  density = 1.86 * g / cm3;
  G4Material* Coating_FEMC = new G4Material("Coating_FEMC", density, ncomponents = 2);
  Coating_FEMC->AddMaterial(G4Material::GetMaterial("Epoxy"), fractionmass = 0.80);
  Coating_FEMC->AddMaterial(material_TiO2, fractionmass = 0.20);

  return Coating_FEMC;
}

//_____________________________________________________________________________
void PHG4LFHcalDetector::SurfaceTable(G4LogicalVolume* vol)
{
  G4OpticalSurface* surface = new G4OpticalSurface("ScintWrapB1");

  new G4LogicalSkinSurface("CrystalSurfaceL", vol, surface);

  surface->SetType(dielectric_metal);
  surface->SetFinish(polished);
  surface->SetModel(glisur);

  //crystal optical surface

  //surface material
  const G4int ntab = 2;
  G4double opt_en[] = {1.551 * eV, 3.545 * eV};  // 350 - 800 nm
  G4double reflectivity[] = {0.9, 0.9};
  G4double efficiency[] = {0.99, 0.99};
  G4MaterialPropertiesTable* surfmat = new G4MaterialPropertiesTable();
  surfmat->AddProperty("REFLECTIVITY", opt_en, reflectivity, ntab);
  surfmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  surface->SetMaterialPropertiesTable(surfmat);
  //csurf->DumpInfo();

}  //SurfaceTable
//_______________________________________________________________________
G4Material* PHG4LFHcalDetector::GetWLSFiberMaterial()
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardEcalDetector: Making WLSFiberFEMC PMMA material..." << std::endl;
  }

  G4double density;
  G4int ncomponents;

  G4Material* material_WLSFiberFEMC = new G4Material("WLSFiberFEMC", density = 1.18 * g / cm3, ncomponents = 3);
  material_WLSFiberFEMC->AddElement(G4Element::GetElement("C"), 5);
  material_WLSFiberFEMC->AddElement(G4Element::GetElement("H"), 8);
  material_WLSFiberFEMC->AddElement(G4Element::GetElement("O"), 2);
  if (m_doLightProp)
  {
    const G4int nEntries = 50;
    G4double photonEnergy[nEntries] =
        {2.00 * eV, 2.03 * eV, 2.06 * eV, 2.09 * eV, 2.12 * eV,
         2.15 * eV, 2.18 * eV, 2.21 * eV, 2.24 * eV, 2.27 * eV,
         2.30 * eV, 2.33 * eV, 2.36 * eV, 2.39 * eV, 2.42 * eV,
         2.45 * eV, 2.48 * eV, 2.51 * eV, 2.54 * eV, 2.57 * eV,
         2.60 * eV, 2.63 * eV, 2.66 * eV, 2.69 * eV, 2.72 * eV,
         2.75 * eV, 2.78 * eV, 2.81 * eV, 2.84 * eV, 2.87 * eV,
         2.90 * eV, 2.93 * eV, 2.96 * eV, 2.99 * eV, 3.02 * eV,
         3.05 * eV, 3.08 * eV, 3.11 * eV, 3.14 * eV, 3.17 * eV,
         3.20 * eV, 3.23 * eV, 3.26 * eV, 3.29 * eV, 3.32 * eV,
         3.35 * eV, 3.38 * eV, 3.41 * eV, 3.44 * eV, 3.47 * eV};
    G4double refractiveIndexWLSfiber[nEntries] =
        {1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
         1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
         1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
         1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
         1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60};

    G4double absWLSfiber[nEntries] =
        {5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
         5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
         5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 1.10 * m,
         1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1. * mm, 1. * mm, 1. * mm, 1. * mm,
         1. * mm, 1. * mm, 1. * mm, 1. * mm, 1. * mm, 1. * mm, 1. * mm, 1. * mm, 1. * mm, 1. * mm};

    G4double emissionFib[nEntries] =
        {0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
         3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
         12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
         15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
         0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
    // Add entries into properties table
    G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
    mptWLSfiber->AddProperty("RINDEX", photonEnergy, refractiveIndexWLSfiber, nEntries);
    mptWLSfiber->AddProperty("WLSABSLENGTH", photonEnergy, absWLSfiber, nEntries);
    mptWLSfiber->AddProperty("WLSCOMPONENT", photonEnergy, emissionFib, nEntries);
    mptWLSfiber->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
    material_WLSFiberFEMC->SetMaterialPropertiesTable(mptWLSfiber);
  }
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardEcalDetector:  Making WLSFiberFEMC material done." << std::endl;
  }

  return material_WLSFiberFEMC;
}

//_____________________________________________________________________________
void PHG4LFHcalDetector::MakeBoundary(G4VPhysicalVolume* crystal, G4VPhysicalVolume* opdet, bool isFiber)
{
  //optical boundary between the crystal and optical photons detector

  G4OpticalSurface* surf = new G4OpticalSurface("OpDetS");
  surf->SetType(dielectric_dielectric);  // photons go to the detector, must have rindex defined
  // surf->SetType(dielectric_metal); // photon is absorbed when reaching the detector, no material rindex required
  //surf->SetFinish(ground);
  surf->SetFinish(polished);
  //surf->SetModel(unified);
  surf->SetModel(glisur);

  new G4LogicalBorderSurface("OpDetB", crystal, opdet, surf);

  const G4int ntab = 2;
  G4double opt_en[] = {1.551 * eV, 3.545 * eV};  // 350 - 800 nm
  //G4double reflectivity[] = {0., 0.};
  G4double reflectivityFiber[] = {0.0, 0.0};
  //G4double reflectivity[] = {0.9, 0.9};
  //G4double reflectivity[] = {1., 1.};
  G4double efficiency[] = {1., 1.};
  G4double RefractiveIndexBoundary[] = {1.0, 1.0};

  G4MaterialPropertiesTable* surfmat = new G4MaterialPropertiesTable();
  surfmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  surfmat->AddProperty("RINDEX", opt_en, RefractiveIndexBoundary, ntab);
  surfmat->AddProperty("REFLECTIVITY", opt_en, reflectivityFiber, ntab);
  surf->SetMaterialPropertiesTable(surfmat);

}  //MakeBoundary

//_____________________________________________________________________________
void PHG4LFHcalDetector::MakeBoundary_Fiber_Scint(G4VPhysicalVolume* fiber, G4VPhysicalVolume* scinti)
{
  //optical boundary between the fiber and scintillator

  G4OpticalSurface* ScintToFiberS = new G4OpticalSurface("ScintToFiberS");
  ScintToFiberS->SetType(dielectric_dielectric);  // photons go to the detector, must have rindex defined
  // ScintToFiberS->SetType(dielectric_metal); // photon is absorbed when reaching the detector, no material rindex required
  //ScintToFiberS->SetFinish(ground);
  ScintToFiberS->SetFinish(groundair);
  // ScintToFiberS->SetFinish(polished);
  //ScintToFiberS->SetModel(unified);
  ScintToFiberS->SetModel(glisur);

  new G4LogicalBorderSurface("ScintToFiberB", scinti, fiber, ScintToFiberS);

  const G4int ntab = 2;
  G4double opt_en[] = {1.551 * eV, 3.545 * eV};  // 350 - 800 nm
  //G4double reflectivity[] = {0., 0.};
  G4double reflectivity[] = {0.01, 0.01};
  //G4double reflectivity[] = {0.9, 0.9};
  //G4double reflectivity[] = {1., 1.};
  G4double efficiency[] = {1., 1.};
  G4double RefractiveIndexBoundary[] = {1.4, 1.4};

  G4MaterialPropertiesTable* ScintToFiberSmat = new G4MaterialPropertiesTable();
  ScintToFiberSmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  ScintToFiberSmat->AddProperty("RINDEX", opt_en, RefractiveIndexBoundary, ntab);
  ScintToFiberSmat->AddProperty("REFLECTIVITY", opt_en, reflectivity, ntab);
  ScintToFiberS->SetMaterialPropertiesTable(ScintToFiberSmat);

  //optical boundary between the fiber and scintillator

  G4OpticalSurface* FiberToScintSurface = new G4OpticalSurface("ScintToFiberS");
  // FiberToScintSurface->SetType(dielectric_dielectric); // photons go to the detector, must have rindex defined
  FiberToScintSurface->SetType(dielectric_metal);  // photon is absorbed when reaching the detector, no material rindex required
  //FiberToScintSurface->SetFinish(ground);
  // FiberToScintSurface->SetFinish(groundair);
  FiberToScintSurface->SetFinish(polished);
  //FiberToScintSurface->SetModel(unified);
  FiberToScintSurface->SetModel(glisur);

  new G4LogicalBorderSurface("ScintToFiberB", scinti, fiber, FiberToScintSurface);

  const G4int ntab2 = 2;
  G4double opt_en2[] = {1.551 * eV, 3.545 * eV};  // 350 - 800 nm
  //G4double reflectivity[] = {0., 0.};
  // G4double reflectivity2[] = {0.01, 0.01};
  //G4double reflectivity[] = {0.9, 0.9};
  G4double reflectivity2[] = {1., 1.};
  G4double efficiency2[] = {1., 1.};
  G4double RefractiveIndexBoundary2[] = {1.6, 1.6};

  G4MaterialPropertiesTable* FiberToScintSurfacemat = new G4MaterialPropertiesTable();
  FiberToScintSurfacemat->AddProperty("EFFICIENCY", opt_en2, efficiency2, ntab2);
  FiberToScintSurfacemat->AddProperty("RINDEX", opt_en2, RefractiveIndexBoundary2, ntab2);
  FiberToScintSurfacemat->AddProperty("REFLECTIVITY", opt_en2, reflectivity2, ntab2);
  FiberToScintSurface->SetMaterialPropertiesTable(FiberToScintSurfacemat);

}  //MakeBoundary

int PHG4LFHcalDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  ifstream istream_mapping;
  istream_mapping.open(m_Params->get_string_param("mapping_file"));
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in PHG4LFHcalDetector: Failed to open mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
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
        std::cout << "PHG4LFHcalDetector: SKIPPING line in mapping file: " << line_mapping << std::endl;
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
        cout << "ERROR in PHG4LFHcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      ostringstream towername;
      towername.str("");
      towername << m_TowerLogicNamePrefix << "_j_" << idx_j << "_k_" << idx_k;

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
      m_TowerPostionMap.insert(make_pair(towername.str(), tower_new));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      string parname;
      G4double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        cout << "ERROR in PHG4LFHcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
        gSystem->Exit(1);
      }

      m_GlobalParameterMap.insert(make_pair(parname, parval));
    }
  }

  /* Update member variables for global parameters based on parsed parameter file */
  std::map<string, G4double>::iterator parit;

  parit = m_GlobalParameterMap.find("Gtower_dx");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("tower_dx", parit->second);  // in cm
  }

  parit = m_GlobalParameterMap.find("Gtower_dy");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("tower_dy", parit->second);  // in cm
  }

  parit = m_GlobalParameterMap.find("Gtower_dz");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("tower_dz", parit->second);  // in cm
  }

  parit = m_GlobalParameterMap.find("Gr1_inner");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rMin1", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gr1_outer");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rMax1", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gr2_inner");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rMin2", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gr2_outer");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rMax2", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gdz");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("dz", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gx0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("place_x", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gy0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("place_y", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gz0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("place_z", parit->second);
  }

  parit = m_GlobalParameterMap.find("Grot_x");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rot_x", parit->second * rad / deg);
  }

  parit = m_GlobalParameterMap.find("Grot_y");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rot_y", parit->second * rad / deg);
  }

  parit = m_GlobalParameterMap.find("Grot_z");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rot_z", parit->second * rad / deg);
  }

  parit = m_GlobalParameterMap.find("thickness_absorber");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("thickness_absorber", parit->second);

  parit = m_GlobalParameterMap.find("thickness_scintillator");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("thickness_scintillator", parit->second);

  parit = m_GlobalParameterMap.find("nlayerspertowerseg");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_int_param("nlayerspertowerseg", (int) parit->second);

  parit = m_GlobalParameterMap.find("xoffset");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("xoffset", parit->second);

  parit = m_GlobalParameterMap.find("yoffset");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("yoffset", parit->second);

  parit = m_GlobalParameterMap.find("frame_width");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("frame_width", parit->second);

  parit = m_GlobalParameterMap.find("width_coating");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("width_coating", parit->second);

  parit = m_GlobalParameterMap.find("usetailcatcher");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_int_param("usetailcatcher", parit->second);

  parit = m_GlobalParameterMap.find("embed_fiber");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_int_param("embed_fiber", parit->second);

  //! TODO make this better!
  if (m_Params->get_int_param("usetailcatcher"))
  {
    m_Params->set_double_param("zdepthcatcheroffset", m_Params->get_double_param("place_z") + (m_Params->get_double_param("tower_dz") / 2) - (10 * (m_Params->get_double_param("thickness_scintillator") + m_Params->get_double_param("thickness_absorber"))));
    m_Params->set_int_param("nLayerOffsetTailcatcher", (int) ((m_Params->get_double_param("tower_dz") - (10 * (m_Params->get_double_param("thickness_scintillator") + m_Params->get_double_param("thickness_absorber")))) / (m_Params->get_double_param("thickness_scintillator") + m_Params->get_double_param("thickness_absorber"))));
  }
  if (Verbosity() > 1)
  {
    std::cout << "PHG4 detector LFHCal - Absorber: " << m_Params->get_double_param("thickness_absorber") << " cm\t Scintilator: " << m_Params->get_double_param("thickness_scintillator") << " cm\t layers per segment: " << m_Params->get_int_param("nlayerspertowerseg") << std::endl;
  }

  return 0;
}
