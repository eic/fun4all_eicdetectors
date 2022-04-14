#include "PHG4ForwardEcalDetector.h"

#include "PHG4ForwardEcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4gdml/PHG4GDMLConfig.hh>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalSkinSurface.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MaterialPropertiesTable.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PVReplica.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4Polyhedra.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <TSystem.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________________
PHG4ForwardEcalDetector::PHG4ForwardEcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4ForwardEcalDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_GdmlConfig(PHG4GDMLUtility::GetOrMakeConfigNode(Node))
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_doLightProp(false)
{
  for (int i = 0; i < 3; i++)
  {
    m_TowerDx[i] = 30 * mm;
    m_TowerDy[i] = 30 * mm;
    m_TowerDz[i] = 170.0 * mm;
  }
  for (int i = 3; i < 7; i++)
  {
    m_TowerDx[i] = NAN;
    m_TowerDy[i] = NAN;
    m_TowerDz[i] = NAN;
  }
  m_RMin[0] = 110 * mm;
  m_RMax[0] = 2250 * mm;
  m_RMin[1] = 120 * mm;
  m_RMax[1] = 2460 * mm;
  m_Params->set_double_param("xoffset", 0.);
  m_Params->set_double_param("yoffset", 0.);

  assert(m_GdmlConfig);
}

//_______________________________________________________________________
int PHG4ForwardEcalDetector::IsInForwardEcal(G4VPhysicalVolume* volume) const
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
void PHG4ForwardEcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardEcalDetector: Begin Construction" << std::endl;
  }

  /* Read parameters for detector construction and mappign from file */
  ParseParametersFromTable();

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));

  G4double tower_readout_dz = m_Params->get_double_param("tower_readout_dz") * cm;

  G4VSolid* beampipe_cutout = new G4Cons("FEMC_beampipe_cutout",
                                         0, m_RMin[0],
                                         0, m_RMin[1],
                                         (m_dZ + tower_readout_dz),
                                         0, 2 * M_PI);
  G4VSolid* ecal_envelope_solid = new G4Cons("FEMC_envelope_solid_cutout",
                                             0, m_RMax[0],
                                             0, m_RMax[1],
                                             (m_dZ + tower_readout_dz) / 2.0,
                                             0, 2 * M_PI);
  ecal_envelope_solid = new G4SubtractionSolid(G4String("hFEMC_envelope_solid"), ecal_envelope_solid, beampipe_cutout, 0, G4ThreeVector(m_Params->get_double_param("xoffset") * cm, m_Params->get_double_param("yoffset") * cm, 0.));

  G4LogicalVolume* ecal_envelope_log = new G4LogicalVolume(ecal_envelope_solid, WorldMaterial, "hFEMC_envelope", 0, 0, 0);

  /* Define visualization attributes for envelope cone */
  GetDisplayAction()->AddVolume(ecal_envelope_log, "Envelope");

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix ecal_rotm;
  ecal_rotm.rotateX(m_XRot);
  ecal_rotm.rotateY(m_YRot);
  ecal_rotm.rotateZ(m_ZRot);

  /* Place envelope cone in simulation */
  std::string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  new G4PVPlacement(G4Transform3D(ecal_rotm, G4ThreeVector(m_PlaceX, m_PlaceY, m_PlaceZ - tower_readout_dz / 2)),
                    ecal_envelope_log, name_envelope, logicWorld, 0, false, OverlapCheck());

  /* Construct single calorimeter towers */
  G4LogicalVolume* singletower[7] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  typedef std::map<std::string, towerposition>::iterator it_type;
  for (it_type iterator = m_TowerPositionMap.begin(); iterator != m_TowerPositionMap.end(); ++iterator)
  {
    for (int i = 0; i < 7; i++)
    {
      if (iterator->second.type == i && singletower[i] == nullptr)
      {
        singletower[i] = ConstructTower(i);
      }
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << singletower << std::endl;
  }
  /* Place calorimeter towers within envelope */
  PlaceTower(ecal_envelope_log, singletower);

  return;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4ForwardEcalDetector::ConstructTower(int type)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardEcalDetector: Build logical volume for single tower, type = " << type << std::endl;
  }
  /* create logical volume for single tower */
  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));

  /* create geometry volumes for scintillator and absorber plates to place inside single_tower */
  // PHENIX EMCal JGL 3/27/2016
  G4int nlayers = 66;
  G4double thickness_layer = m_TowerDz[2] / (float) nlayers;
  // update layer thickness with https://doi.org/10.1016/S0168-9002(02)01954-X
  G4double thickness_cell = 5.6 * mm;
  G4double thickness_absorber = 1.55 * mm;     // 1.55mm absorber
  G4double thickness_scintillator = 4.0 * mm;  // 4mm scintillator
  // notched in TiO2 coating in scintillator plate
  G4double width_coating = m_Params->get_double_param("width_coating") * cm;  // 4mm scintillator
  G4Material* material_scintillator = GetScintillatorMaterial();              //G4Material::GetMaterial("G4_POLYSTYRENE");
  G4Material* material_absorber = GetDetectorMaterial("G4_Pb");
  G4int nFibers = m_Params->get_int_param("nFibers");
  G4double fiber_diam = m_Params->get_double_param("fiber_diam") * cm;  // 4mm scintillator
  G4double fiber_extra_length = 0.0;
  // additional fiber length at end of calorimeter
  if (fiber_diam > 0) fiber_extra_length = 0.5 * cm;
  // width of plate in front of FEMC for clamping all layers
  G4double clamp_plate_width = m_Params->get_double_param("clamp_plate_width") * cm;
  // depth of readout (4cm)
  G4double tower_readout_dz = m_Params->get_double_param("tower_readout_dz") * cm;

  G4VSolid* single_tower_solid = new G4Box("single_tower_solid2",
                                           m_TowerDx[2] / 2.0,
                                           m_TowerDy[2] / 2.0,
                                           (m_TowerDz[2] + tower_readout_dz) / 2.0);
  G4VSolid* single_tower_solid_replica = new G4Box("single_tower_solid_replica",
                                                   m_TowerDx[2] / 2.0,
                                                   m_TowerDy[2] / 2.0,
                                                   m_TowerDz[2] / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            WorldMaterial,
                                                            "single_tower_logic2",
                                                            0, 0, 0);

  GetDisplayAction()->AddVolume(single_tower_logic, "SingleTower");

  if (Verbosity())
  {
    std::cout << " m_TowerDz[2] = " << m_TowerDz[2] << " thickness_layer = " << thickness_layer << " thickness_cell = " << thickness_cell << std::endl;
  }

  if (thickness_layer <= thickness_cell)
  {
    std::cout << __PRETTY_FUNCTION__
              << "Tower size z (m_TowerDz[2) from database is too thin. "
              << "It does not fit the layer structure as described in https://doi.org/10.1016/S0168-9002(02)01954-X !" << std::endl
              << "Abort" << std::endl;
    std::cout << " m_TowerDz[2] = " << m_TowerDz[2] << " i.e. nlayers " << nlayers << " * thickness_layer " << thickness_layer << " <= thickness_cell " << thickness_cell << std::endl;
    exit(1);
  }

  //**********************************************************************************************
  /* create logical and geometry volumes for minitower read-out unit */
  //**********************************************************************************************
  G4VSolid* miniblock_solid = new G4Box("miniblock_solid",
                                        m_TowerDx[2] / 2.0,
                                        m_TowerDy[2] / 2.0,
                                        thickness_cell / 2.0);
  G4LogicalVolume* miniblock_logic = new G4LogicalVolume(miniblock_solid,
                                                         WorldMaterial,
                                                         "miniblock_logic",
                                                         0, 0, 0);
  GetDisplayAction()->AddVolume(miniblock_logic, "miniblock");
  //**********************************************************************************************
  /* create logical & geometry volumes for scintillator and absorber plates to place inside mini read-out unit */
  //**********************************************************************************************
  G4VSolid* solid_absorber = new G4Box("single_plate_absorber_solid2",
                                       m_TowerDx[2] / 2.0,
                                       m_TowerDy[2] / 2.0,
                                       thickness_absorber / 2.0);

  G4VSolid* solid_scintillator = new G4Box("single_plate_scintillator2",
                                           (m_TowerDx[2]) / 2.0,
                                           (m_TowerDy[2]) / 2.0,
                                           //  (m_TowerDx[2] - 2*width_coating) / 2.0,
                                           //  (m_TowerDy[2] - 2*width_coating) / 2.0,
                                           thickness_scintillator / 2.0);

  if (clamp_plate_width > 0)
  {
    G4VSolid* solid_clamp1 = new G4Box("single_plate_clamp_solid1",
                                       m_TowerDx[2] / 2.0,
                                       m_TowerDy[2] / 2.0,
                                       clamp_plate_width / 2.0);
    G4LogicalVolume* logic_clampplate = new G4LogicalVolume(solid_clamp1,
                                                            GetDetectorMaterial("G4_Fe"),
                                                            "logic_clampplate",
                                                            0, 0, 0);
    m_AbsorberLogicalVolSet.insert(logic_clampplate);
    GetDisplayAction()->AddVolume(logic_clampplate, "Clamp");
    std::string name_clamp = m_TowerLogicNamePrefix + "_single_plate_clamp";

    new G4PVPlacement(0, G4ThreeVector(0, 0, (m_dZ + tower_readout_dz) / 2.0 - clamp_plate_width / 2.0),
                      logic_clampplate,
                      name_clamp,
                      single_tower_logic,
                      0, 0, OverlapCheck());
  }
  if (nFibers > 0 && nFibers == 5 && fiber_diam > 0)
  {
    G4VSolid* cutoutfiber_solid = new G4Tubs("cutoutfiber_solid",
                                             0.0, 1.01 * fiber_diam / 2.0, m_TowerDz[2], 0.0, 2 * M_PI);

    single_tower_solid_replica = new G4SubtractionSolid(G4String("single_tower_solid_replica_cu1"), single_tower_solid_replica, cutoutfiber_solid, 0, G4ThreeVector(0, 0, 0.));
    single_tower_solid_replica = new G4SubtractionSolid(G4String("single_tower_solid_replica_cu2"), single_tower_solid_replica, cutoutfiber_solid, 0, G4ThreeVector(-m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, 0.));
    single_tower_solid_replica = new G4SubtractionSolid(G4String("single_tower_solid_replica_cu3"), single_tower_solid_replica, cutoutfiber_solid, 0, G4ThreeVector(m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, 0.));
    single_tower_solid_replica = new G4SubtractionSolid(G4String("single_tower_solid_replica_cu4"), single_tower_solid_replica, cutoutfiber_solid, 0, G4ThreeVector(m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, 0.));
    single_tower_solid_replica = new G4SubtractionSolid(G4String("single_tower_solid_replica_cu5"), single_tower_solid_replica, cutoutfiber_solid, 0, G4ThreeVector(-m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, 0.));

    solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_cu1"), solid_absorber, cutoutfiber_solid, 0, G4ThreeVector(0, 0, 0.));
    solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_cu2"), solid_absorber, cutoutfiber_solid, 0, G4ThreeVector(-m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, 0.));
    solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_cu3"), solid_absorber, cutoutfiber_solid, 0, G4ThreeVector(m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, 0.));
    solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_cu4"), solid_absorber, cutoutfiber_solid, 0, G4ThreeVector(m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, 0.));
    solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_cu5"), solid_absorber, cutoutfiber_solid, 0, G4ThreeVector(-m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, 0.));

    solid_scintillator = new G4SubtractionSolid(G4String("solid_scintillator_cu1"), solid_scintillator, cutoutfiber_solid, 0, G4ThreeVector(0, 0, 0.));
    solid_scintillator = new G4SubtractionSolid(G4String("solid_scintillator_cu2"), solid_scintillator, cutoutfiber_solid, 0, G4ThreeVector(-m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, 0.));
    solid_scintillator = new G4SubtractionSolid(G4String("solid_scintillator_cu3"), solid_scintillator, cutoutfiber_solid, 0, G4ThreeVector(m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, 0.));
    solid_scintillator = new G4SubtractionSolid(G4String("solid_scintillator_cu4"), solid_scintillator, cutoutfiber_solid, 0, G4ThreeVector(m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, 0.));
    solid_scintillator = new G4SubtractionSolid(G4String("solid_scintillator_cu5"), solid_scintillator, cutoutfiber_solid, 0, G4ThreeVector(-m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, 0.));

    if (clamp_plate_width > 0)
    {
      G4VSolid* solid_clamp2 = new G4Box("single_plate_clamp_solid2",
                                         m_TowerDx[2] / 2.0,
                                         m_TowerDy[2] / 2.0,
                                         clamp_plate_width / 2.0);
      solid_clamp2 = new G4SubtractionSolid(G4String("solid_clamp2_cu1"), solid_clamp2, cutoutfiber_solid, 0, G4ThreeVector(0, 0, 0.));
      solid_clamp2 = new G4SubtractionSolid(G4String("solid_clamp2_cu2"), solid_clamp2, cutoutfiber_solid, 0, G4ThreeVector(-m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, 0.));
      solid_clamp2 = new G4SubtractionSolid(G4String("solid_clamp2_cu3"), solid_clamp2, cutoutfiber_solid, 0, G4ThreeVector(m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, 0.));
      solid_clamp2 = new G4SubtractionSolid(G4String("solid_clamp2_cu4"), solid_clamp2, cutoutfiber_solid, 0, G4ThreeVector(m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, 0.));
      solid_clamp2 = new G4SubtractionSolid(G4String("solid_clamp2_cu5"), solid_clamp2, cutoutfiber_solid, 0, G4ThreeVector(-m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, 0.));
      G4LogicalVolume* logic_clampplate2 = new G4LogicalVolume(solid_clamp2,
                                                               GetDetectorMaterial("G4_Fe"),
                                                               "logic_clampplate2",
                                                               0, 0, 0);
      m_AbsorberLogicalVolSet.insert(logic_clampplate2);
      GetDisplayAction()->AddVolume(logic_clampplate2, "Clamp");
      std::string name_clamp = m_TowerLogicNamePrefix + "_single_plate_clamp2";

      new G4PVPlacement(0, G4ThreeVector(0, 0, (tower_readout_dz) / 2.0 - clamp_plate_width - m_TowerDz[2] / 2.0 - clamp_plate_width / 2.0),
                        logic_clampplate2,
                        name_clamp,
                        single_tower_logic,
                        0, 0, OverlapCheck());
    }
  }

  if (width_coating > 0)
  {
    G4double depthCoating = 0.93 * thickness_scintillator;
    G4double zPlaneCoating[2] = {0, depthCoating};
    G4double rInnerCoating[2] = {(m_TowerDx[2] - 2 * width_coating) / 2, (m_TowerDx[2] - 2 * width_coating) / 2.0};
    G4double rOuterCoating[2] = {m_TowerDx[2] / 2.0, m_TowerDx[2] / 2.0};
    G4VSolid* solid_coating = new G4Polyhedra("solid_coating",
                                              0, 2 * M_PI,
                                              4, 2,
                                              zPlaneCoating, rInnerCoating, rOuterCoating);
    G4LogicalVolume* logic_coating = new G4LogicalVolume(solid_coating,
                                                         GetCoatingMaterial(),
                                                         "logic_coating",
                                                         0, 0, 0);
    m_AbsorberLogicalVolSet.insert(logic_coating);
    GetDisplayAction()->AddVolume(logic_coating, "Coating");
    std::string name_coating = m_TowerLogicNamePrefix + "_single_plate_coating";

    G4RotationMatrix* rotCoating = new G4RotationMatrix();
    rotCoating->rotateZ(M_PI / 4);
    new G4PVPlacement(rotCoating, G4ThreeVector(0, 0, thickness_cell / 2.0 - depthCoating),
                      logic_coating,
                      name_coating,
                      miniblock_logic,
                      0, 0, OverlapCheck());
    solid_scintillator = new G4SubtractionSolid(G4String("solid_scintillator_cuCoating"), solid_scintillator, solid_coating, rotCoating, G4ThreeVector(0, 0, thickness_scintillator / 2.0 - depthCoating));
  }
  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "single_plate_absorber_logic2",
                                                        0, 0, 0);
  m_AbsorberLogicalVolSet.insert(logic_absorber);
  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                     material_scintillator,
                                                     "hEcal_scintillator_plate_logic2",
                                                     0, 0, 0);
  m_ScintiLogicalVolSet.insert(logic_scint);
  if (m_doLightProp)
  {
    SurfaceTable(logic_scint);
  }

  GetDisplayAction()->AddVolume(logic_absorber, "Absorber");
  GetDisplayAction()->AddVolume(logic_scint, "Scintillator");

  std::string name_absorber = m_TowerLogicNamePrefix + "_single_plate_absorber2";
  std::string name_scintillator = m_TowerLogicNamePrefix + "_single_plate_scintillator2";

  new G4PVPlacement(0, G4ThreeVector(0, 0, -thickness_cell / 2.0 + thickness_absorber / 2.0),
                    logic_absorber,
                    name_absorber,
                    miniblock_logic,
                    0, 0, OverlapCheck());

  new G4PVPlacement(0, G4ThreeVector(0, 0, thickness_cell / 2.0 - thickness_scintillator / 2.0),
                    logic_scint,
                    name_scintillator,
                    miniblock_logic,
                    0, 0, OverlapCheck());

  G4LogicalVolume* single_tower_replica_logic = new G4LogicalVolume(single_tower_solid_replica,
                                                                    WorldMaterial,
                                                                    "single_tower_replica_logic",
                                                                    0, 0, 0);
  /* create replica within tower */
  std::string name_tower = m_TowerLogicNamePrefix;
  new G4PVReplica(name_tower, miniblock_logic, single_tower_replica_logic,
                  kZAxis, nlayers, thickness_layer, 0);

  GetDisplayAction()->AddVolume(single_tower_replica_logic, "SingleTower");

  new G4PVPlacement(0, G4ThreeVector(0, 0, (tower_readout_dz) / 2.0 - clamp_plate_width),
                    single_tower_replica_logic,
                    "replicated_layers_placed",
                    single_tower_logic,
                    0, 0, OverlapCheck());

  // place array of fibers inside absorber
  if (nFibers > 0 && nFibers == 5 && fiber_diam > 0)
  {
    std::string fiberName = "single_fiber_scintillator_solid" + std::to_string(type);
    G4VSolid* single_scintillator_solid = new G4Tubs(fiberName,
                                                     0.0, fiber_diam / 2.0, (m_TowerDz[2] + fiber_extra_length) / 2, 0.0, 2 * M_PI);

    /* create logical volumes for scintillator and absorber plates to place inside single_tower */
    G4Material* material_WLSFiber = GetWLSFiberFEMCMaterial();

    std::string fiberLogicName = "hEcal_scintillator_fiber_logic" + std::to_string(type);
    G4LogicalVolume* single_scintillator_logic = new G4LogicalVolume(single_scintillator_solid,
                                                                     material_WLSFiber,
                                                                     fiberLogicName,
                                                                     0, 0, 0);
    m_ScintiLogicalVolSet.insert(single_scintillator_logic);
    GetDisplayAction()->AddVolume(single_scintillator_logic, "Fiber");

    std::string name_scintillator = m_TowerLogicNamePrefix + "_single_fiber_scintillator" + std::to_string(type);

    new G4PVPlacement(0, G4ThreeVector(0, 0, -fiber_extra_length / 2 + tower_readout_dz / 2 - clamp_plate_width),
                      single_scintillator_logic,
                      name_scintillator + "_center",
                      single_tower_logic,
                      0, 0, OverlapCheck());

    new G4PVPlacement(0, G4ThreeVector(-m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, -fiber_extra_length / 2 + tower_readout_dz / 2 - clamp_plate_width),
                      single_scintillator_logic,
                      name_scintillator + "_tl",
                      single_tower_logic,
                      0, 0, OverlapCheck());

    new G4PVPlacement(0, G4ThreeVector(m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, -fiber_extra_length / 2 + tower_readout_dz / 2 - clamp_plate_width),
                      single_scintillator_logic,
                      name_scintillator + "_bl",
                      single_tower_logic,
                      0, 0, OverlapCheck());

    new G4PVPlacement(0, G4ThreeVector(m_TowerDx[2] / 4.0, m_TowerDy[2] / 4.0, -fiber_extra_length / 2 + tower_readout_dz / 2 - clamp_plate_width),
                      single_scintillator_logic,
                      name_scintillator + "_tr",
                      single_tower_logic,
                      0, 0, OverlapCheck());

    new G4PVPlacement(0, G4ThreeVector(-m_TowerDx[2] / 4.0, -m_TowerDy[2] / 4.0, -fiber_extra_length / 2 + tower_readout_dz / 2 - clamp_plate_width),
                      single_scintillator_logic,
                      name_scintillator + "_br",
                      single_tower_logic,
                      0, 0, OverlapCheck());
  }

  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardEcalDetector: Building logical volume for single tower done." << std::endl;
  }

  return single_tower_logic;
}

int PHG4ForwardEcalDetector::PlaceTower(G4LogicalVolume* ecalenvelope, G4LogicalVolume* singletowerIn[7])
{
  /* Loop over all tower positions in vector and place tower */
  for (std::map<std::string, towerposition>::iterator iterator = m_TowerPositionMap.begin(); iterator != m_TowerPositionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4ForwardEcalDetector: Place tower " << iterator->first
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k
                << " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << std::endl;
    }

    assert(iterator->second.type >= 0 && iterator->second.type <= 6);
    G4LogicalVolume* singletower = singletowerIn[iterator->second.type];
    int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;

    G4PVPlacement* tower_placement =
        new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                          singletower,
                          iterator->first,
                          ecalenvelope,
                          0, copyno, OverlapCheck());

    m_GdmlConfig->exclude_physical_vol(tower_placement);
  }

  return 0;
}

//_______________________________________________________________________
G4Material* PHG4ForwardEcalDetector::GetScintillatorMaterial()
{
  G4double density;
  G4int ncomponents;
  G4Material* material_ScintFEMC = new G4Material("PolystyreneFEMC", density = 1.03 * g / cm3, ncomponents = 2);
  material_ScintFEMC->AddElement(G4Element::GetElement("C"), 8);
  material_ScintFEMC->AddElement(G4Element::GetElement("H"), 8);

  if (m_doLightProp)
  {
    const G4int ntab = 4;

    G4double wls_Energy[] = {2.00 * eV, 2.87 * eV, 2.90 * eV,
                             3.47 * eV};

    G4double rIndexPstyrene[] = {1.5, 1.5, 1.5, 1.5};
    G4double absorption1[] = {2. * cm, 2. * cm, 2. * cm, 2. * cm};
    G4double scintilFast[] = {0.0, 0.0, 1.0, 1.0};
    G4MaterialPropertiesTable* fMPTPStyrene = new G4MaterialPropertiesTable();
    fMPTPStyrene->AddProperty("RINDEX", wls_Energy, rIndexPstyrene, ntab);
    fMPTPStyrene->AddProperty("ABSLENGTH", wls_Energy, absorption1, ntab);
    fMPTPStyrene->AddProperty("SCINTILLATIONCOMPONENT1", wls_Energy, scintilFast, ntab);
    fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
    fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fMPTPStyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT", 10. * ns);

    material_ScintFEMC->SetMaterialPropertiesTable(fMPTPStyrene);
  }
  // Set the Birks Constant for the Polystyrene scintillator
  material_ScintFEMC->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  return material_ScintFEMC;
}

//_______________________________________________________________________
G4Material* PHG4ForwardEcalDetector::GetCoatingMaterial()
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
void PHG4ForwardEcalDetector::SurfaceTable(G4LogicalVolume* vol)
{
  G4OpticalSurface* surface = new G4OpticalSurface("ScintWrapB1");

  new G4LogicalSkinSurface("CrystalSurfaceL", vol, surface);

  surface->SetType(dielectric_metal);
  surface->SetFinish(polished);
  surface->SetModel(glisur);

  //crystal optical surface

  //surface material
  // const G4int ntab = 2;
  // G4double opt_en[] = {1.551*eV, 3.545*eV}; // 350 - 800 nm
  // G4double reflectivity[] = {0.8, 0.8};
  // G4double efficiency[] = {0.9, 0.9};
  G4MaterialPropertiesTable* surfmat = new G4MaterialPropertiesTable();
  // surfmat->AddProperty("REFLECTIVITY", opt_en, reflectivity, ntab);
  // surfmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  surface->SetMaterialPropertiesTable(surfmat);
  //csurf->DumpInfo();

}  //SurfaceTable
//_______________________________________________________________________
G4Material* PHG4ForwardEcalDetector::GetWLSFiberFEMCMaterial()
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardEcalDetector: Making WLSFiberFEMC material..." << std::endl;
  }

  G4double density;
  G4int ncomponents;

  G4Material* material_WLSFiberFEMC = new G4Material("WLSFiberFEMC", density = 1.18 * g / cm3, ncomponents = 3);
  material_WLSFiberFEMC->AddElement(G4Element::GetElement("C"), 5);
  material_WLSFiberFEMC->AddElement(G4Element::GetElement("H"), 8);
  material_WLSFiberFEMC->AddElement(G4Element::GetElement("O"), 2);
  if (m_doLightProp)
  {
    const G4int ntab = 4;
    G4double wls_Energy[] = {2.00 * eV, 2.87 * eV, 2.90 * eV,
                             3.47 * eV};

    G4double RefractiveIndexFiber[] = {1.6, 1.6, 1.6, 1.6};
    G4double AbsFiber[] = {9.0 * m, 9.0 * m, 0.1 * mm, 0.1 * mm};
    G4double EmissionFib[] = {1.0, 1.0, 0.0, 0.0};
    // Add entries into properties table
    G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
    mptWLSfiber->AddProperty("RINDEX", wls_Energy, RefractiveIndexFiber, ntab);
    mptWLSfiber->AddProperty("WLSABSLENGTH", wls_Energy, AbsFiber, ntab);
    mptWLSfiber->AddProperty("WLSCOMPONENT", wls_Energy, EmissionFib, ntab);
    mptWLSfiber->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
    material_WLSFiberFEMC->SetMaterialPropertiesTable(mptWLSfiber);
  }
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardEcalDetector:  Making WLSFiberFEMC material done." << std::endl;
  }

  return material_WLSFiberFEMC;
}

int PHG4ForwardEcalDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  std::ifstream istream_mapping;
  istream_mapping.open(m_Params->get_string_param("mapping_file"));
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in PHG4ForwardEcalDetector: Failed to open mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
    gSystem->Exit(1);
  }

  /* loop over lines in file */
  std::string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {
    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != std::string::npos)
    {
      if (Verbosity() > 0)
      {
        std::cout << "PHG4ForwardEcalDetector: SKIPPING line in mapping file: " << line_mapping << std::endl;
      }
      continue;
    }

    std::istringstream iss(line_mapping);

    /* If line starts with keyword Tower, add to tower positions */
    if (line_mapping.find("Tower ") != std::string::npos)
    {
      unsigned idx_j, idx_k, idx_l;
      G4double pos_x, pos_y, pos_z;
      G4double size_x, size_y, size_z;
      G4double rot_x, rot_y, rot_z;
      int type;
      std::string dummys;

      /* read string- break if error */
      if (!(iss >> dummys >> type >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z))
      {
        std::cout << "ERROR in PHG4ForwardEcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      std::ostringstream towername;
      towername << m_TowerLogicNamePrefix << "_t_" << type << "_j_" << idx_j << "_k_" << idx_k;
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
      tower_new.type = type;
      m_TowerPositionMap.insert(std::make_pair(towername.str(), tower_new));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      std::string parname;
      double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        std::cout << "ERROR in PHG4ForwardEcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
        gSystem->Exit(1);
      }

      m_GlobalParameterMap.insert(std::make_pair(parname, parval));
    }
  }
  /* Update member variables for global parameters based on parsed parameter file */
  std::map<std::string, double>::iterator parit;
  std::ostringstream twr;
  for (int i = 0; i < 7; i++)
  {
    twr.str("");
    twr << "Gtower" << i << "_dx";
    parit = m_GlobalParameterMap.find(twr.str());
    m_TowerDx[i] = parit->second * cm;
    twr.str("");
    twr << "Gtower" << i << "_dy";
    parit = m_GlobalParameterMap.find(twr.str());
    m_TowerDy[i] = parit->second * cm;
    twr.str("");
    twr << "Gtower" << i << "_dz";
    parit = m_GlobalParameterMap.find(twr.str());
    m_TowerDz[i] = parit->second * cm;
  }
  std::ostringstream rad;
  for (int i = 0; i < 2; i++)
  {
    int index = i + 1;
    rad.str("");
    rad << "Gr" << index << "_inner";
    parit = m_GlobalParameterMap.find(rad.str());
    if (parit != m_GlobalParameterMap.end())
    {
      m_RMin[i] = parit->second * cm;
    }
    rad.str("");
    rad << "Gr" << index << "_outer";
    parit = m_GlobalParameterMap.find(rad.str());
    if (parit != m_GlobalParameterMap.end())
    {
      m_RMax[i] = parit->second * cm;
    }
  }
  parit = m_GlobalParameterMap.find("Gdz");
  if (parit != m_GlobalParameterMap.end())
  {
    m_dZ = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Gx0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_PlaceX = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Gy0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_PlaceY = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Gz0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_PlaceZ = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Grot_x");
  if (parit != m_GlobalParameterMap.end())
  {
    m_XRot = parit->second;
  }
  parit = m_GlobalParameterMap.find("Grot_y");
  if (parit != m_GlobalParameterMap.end())
  {
    m_YRot = parit->second;
  }
  parit = m_GlobalParameterMap.find("Grot_z");
  if (parit != m_GlobalParameterMap.end())
  {
    m_ZRot = parit->second;
  }
  parit = m_GlobalParameterMap.find("tower_type");
  if (parit != m_GlobalParameterMap.end())
  {
    m_TowerType = parit->second;
  }

  parit = m_GlobalParameterMap.find("xoffset");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("xoffset", parit->second);

  parit = m_GlobalParameterMap.find("yoffset");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("yoffset", parit->second);

  parit = m_GlobalParameterMap.find("width_coating");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("width_coating", parit->second);

  parit = m_GlobalParameterMap.find("clamp_plate_width");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("clamp_plate_width", parit->second);

  parit = m_GlobalParameterMap.find("fiber_diam");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_double_param("fiber_diam", parit->second);

  parit = m_GlobalParameterMap.find("nFibers");
  if (parit != m_GlobalParameterMap.end())
    m_Params->set_int_param("nFibers", parit->second);

  return 0;
}

void PHG4ForwardEcalDetector::SetTowerDimensions(double dx, double dy, double dz, int type)
{
  assert(type >= 0 && type <= 6);
  m_TowerDx[type] = dx;
  m_TowerDy[type] = dy;
  m_TowerDz[type] = dz;
}
