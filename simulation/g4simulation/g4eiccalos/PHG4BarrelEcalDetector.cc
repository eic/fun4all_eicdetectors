#include "PHG4BarrelEcalDetector.h"
#include "PHG4BarrelEcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <TSystem.h>
#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>
#include <phool/recoConsts.h>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4CutTubs.hh>
#include <Geant4/G4DisplacedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MaterialPropertiesTable.hh>  // for G4MaterialProperties...
#include <Geant4/G4MaterialPropertyVector.hh>   // for G4MaterialPropertyVector
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4TwoVector.hh>  // for G4VPhysicalVolume
#include <Geant4/G4GenericTrap.hh>  // for G4VPhysicalVolume

#include <g4gdml/PHG4GDMLConfig.hh>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

using namespace std;

class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________________
PHG4BarrelEcalDetector::PHG4BarrelEcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4BarrelEcalDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_SupportActiveFlag(m_Params->get_int_param("supportactive"))
  , m_TowerLogicNamePrefix("bcalTower")
  , m_SuperDetector("NONE")
  , m_useLeadGlass(false)
{
  gdml_config = PHG4GDMLUtility::GetOrMakeConfigNode(Node);
  assert(gdml_config);
}
//_______________________________________________________________________
int PHG4BarrelEcalDetector::IsInBarrelEcal(G4VPhysicalVolume* volume) const
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

  if (m_SupportActiveFlag)
  {
    if (m_SupportLogicalVolSet.find(mylogvol) != m_SupportLogicalVolSet.end())
    {
      return -2;
    }
  }
  return 0;
}

//_______________________________________________________________________
void PHG4BarrelEcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4BarrelEcalDetector: Begin Construction" << std::endl;
  }

  if (m_Params->get_string_param("mapping_file").empty())
  {
    std::cout << "ERROR in PHG4BarrelEcalDetector: No mapping file specified. Abort detector construction." << std::endl;
    std::cout << "Please run set_string_param(\"mapping_file\", std::string filename ) first." << std::endl;
    gSystem->Exit(1);
  }

  ParseParametersFromTable();

  G4double Radius = m_Params->get_double_param("radius") * cm;
  G4double tower_length = m_Params->get_double_param("tower_length") * cm;
  G4double becal_length = m_Params->get_double_param("becal_length") * cm;
  G4double CenterZ_Shift = m_Params->get_double_param("CenterZ_Shift") * cm;
  G4double cone1_h = m_Params->get_double_param("cone1_h") * cm;
  G4double cone1_dz = m_Params->get_double_param("cone1_dz") * cm;
  G4double cone2_h = m_Params->get_double_param("cone2_h") * cm;
  G4double cone2_dz = m_Params->get_double_param("cone2_dz") * cm;

  silicon_width_half = m_Params->get_double_param("silicon_width_half") * cm;
  kapton_width_half = m_Params->get_double_param("kapton_width_half") * cm;
  SIO2_width_half = m_Params->get_double_param("SIO2_width_half") * cm;
  Carbon_width_half = m_Params->get_double_param("Carbon_width_half") * cm;
  support_length = m_Params->get_double_param("support_length") * cm;

  G4double max_radius = Radius + tower_length + 2 * silicon_width_half + 2 * kapton_width_half + 2 * SIO2_width_half + 2 * Carbon_width_half + support_length + 8 * overlap;

  //std::cout << Radius << "  " << max_radius << "====================================" << std::endl;

  G4double pos_x1 = 0 * cm;
  G4double pos_y1 = 0 * cm;
  G4double pos_z1 = 0 * cm;

  G4Tubs* cylinder_solid1 = new G4Tubs("BCAL_SOLID1",
                                       Radius, max_radius,
                                       becal_length / 2, 0, 2 * M_PI);

  G4Tubs* cylinder_solid2 = new G4Tubs("BCAL_SOLID2",
                                       Radius - 1, max_radius + 1,
                                       abs(CenterZ_Shift), 0, 2 * M_PI);

  G4ThreeVector shift_cs2 = G4ThreeVector(0, 0, becal_length / 2);

  G4VSolid* cylinder_solid3 = new G4SubtractionSolid("BCAL_SOLID3", cylinder_solid1, cylinder_solid2, 0, shift_cs2);

  G4Cons* cone1 = new G4Cons("cone1",
                             Radius - 1, Radius - 1,
                             Radius - 1, Radius + cone1_h,
                             cone1_dz, 0, 2 * M_PI);

  G4ThreeVector shift_cone1 = G4ThreeVector(0, 0, becal_length / 2 - abs(CenterZ_Shift) - cone1_dz);

  G4VSolid* cylinder_solid4 = new G4SubtractionSolid("BCAL_SOLID4", cylinder_solid3, cone1, 0, shift_cone1);

  G4Cons* cone2 = new G4Cons("cone2",
                             Radius - 1, Radius + cone2_h,
                             Radius - 1, Radius - 1,
                             cone2_dz, 0, 2 * M_PI);

  G4ThreeVector shift_cone2 = G4ThreeVector(0, 0, -becal_length / 2 + cone2_dz);

  G4VSolid* cylinder_solid = new G4SubtractionSolid("BCAL_SOLID", cylinder_solid4, cone2, 0, shift_cone2);

  G4Material* cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);

  G4LogicalVolume* cylinder_logic = new G4LogicalVolume(cylinder_solid, cylinder_mat,
                                                        "BCAL_SOLID", 0, 0, 0);

  m_DisplayAction->AddVolume(cylinder_logic, "BCalCylinder");

  std::string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  G4PVPlacement* phys_envelope =
      new G4PVPlacement(0, G4ThreeVector(pos_x1, pos_y1, pos_z1), cylinder_logic, name_envelope,
                        logicWorld, false, 0, OverlapCheck());

  gdml_config->exclude_physical_vol(phys_envelope);
  PlaceTower(cylinder_logic);

  return;
}

int PHG4BarrelEcalDetector::PlaceTower(G4LogicalVolume* sec)
{
  /* Loop over all tower positions in vector and place tower */
  for (std::map<std::string, towerposition>::iterator iterator = m_TowerPostionMap.begin(); iterator != m_TowerPostionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4BarrelEcalDetector: Place tower " << iterator->first
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k << std::endl;
    }

    int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;

    G4LogicalVolume* block_logic = ConstructTower(iterator);
    // m_DisplayAction->AddVolume(block_logic, iterator->first);
    if (iterator->second.idx_k % 2 == 0)
    {
      m_DisplayAction->AddVolume(block_logic, "Invisible");
    }
    else
    {
      m_DisplayAction->AddVolume(block_logic, "Invisible");
    }

    G4RotationMatrix becal_rotm;
    becal_rotm.rotateX(iterator->second.rotx);
    becal_rotm.rotateY(iterator->second.roty);
    becal_rotm.rotateZ(iterator->second.rotz);

    new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(iterator->second.centerx, iterator->second.centery, iterator->second.centerz)),
                      block_logic,
                      G4String(string(iterator->first) + string("_TT")),
                      sec,
                      0, copyno, OverlapCheck());

  }
  return 0;
}

G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructTower(std::map<std::string, towerposition>::iterator iterator)
{
  G4GenericTrap* block_tower = GetTowerTrap(iterator);
  G4GenericTrap* block_glass = GetGlassTrap(iterator,false);
  G4GenericTrap* block_glass_subtract = GetGlassTrap(iterator,true);
  G4VSolid* block_shell = new G4SubtractionSolid(G4String(string(iterator->first) + string("_Envelope")), block_tower, block_glass_subtract, 0, G4ThreeVector(0.,0.,0.));


  G4Material* material_scinti = GetSciGlass();
  assert(material_scinti);
  G4Material* material_shell = GetCarbonFiber();
  assert(material_shell);

  // G4Trap* block_tower = GetTowerTrap(iterator);
  // G4Trap* block_glass = GetGlassTrapSubtract(iterator);
  // G4ThreeVector shift = G4ThreeVector(-th * sin(iterator->second.roty - M_PI_2), 0, th / 2 * abs(cos(iterator->second.roty - M_PI_2)));
  // G4VSolid* block_solid = new G4SubtractionSolid(G4String(string(iterator->first) + string("_Envelope")), block_tower, block_glass, 0, shift);
  // G4Material* material_shell = GetCarbonFiber();
  // assert(material_shell);

  // G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, material_shell,
  //                                                    G4String(string(iterator->first) + string("_Tower")), 0, 0,
  //                                                    nullptr);
  // m_AbsorberLogicalVolSet.insert(block_logic);


  G4LogicalVolume* block_logic = new G4LogicalVolume(block_tower, G4Material::GetMaterial("G4_AIR"),
                                                    G4String(string(iterator->first) + string("_Tower")), 0, 0,
                                                    nullptr);
  
  G4LogicalVolume* glass_logic = new G4LogicalVolume(block_glass, material_scinti,
                                                    G4String(string(iterator->first) + string("_Glass")), 0, 0,
                                                    nullptr);
  
  G4LogicalVolume* shell_logic = new G4LogicalVolume(block_shell, material_shell,
                                                    G4String(string(iterator->first) + string("_Wall")), 0, 0,
                                                    nullptr);
  new G4PVPlacement(0,G4ThreeVector(0, 0, 0),
                  glass_logic,
                  G4String(string(iterator->first) + string("_Glass")),
                  block_logic,
                  0, 0, OverlapCheck());
  new G4PVPlacement(0,G4ThreeVector(0, 0, 0),
                  shell_logic,
                  G4String(string(iterator->first) + string("_Shell")),
                  block_logic,
                  0, 0, OverlapCheck());

  if (iterator->second.idx_k % 2 == 0)
  {
    m_DisplayAction->AddVolume(glass_logic, "Block1");
  }
  else
  {
    m_DisplayAction->AddVolume(glass_logic, "Block2");
  }
  m_DisplayAction->AddVolume(shell_logic, "Carbon");

  m_AbsorberLogicalVolSet.insert(shell_logic);
  m_ScintiLogicalVolSet.insert(glass_logic);
  return block_logic;
}

G4Material* PHG4BarrelEcalDetector::GetCarbonFiber()
{
  static string matname = "CrystalCarbonFiber";
  G4Material* carbonfiber = G4Material::GetMaterial(matname, false);  // false suppresses warning that material does not exist
  if (!carbonfiber)
  {
    G4double density_carbon_fiber = 1.44 * g / cm3;
    carbonfiber = new G4Material(matname, density_carbon_fiber, 1);
    carbonfiber->AddElement(G4Element::GetElement("C"), 1);
  }
  return carbonfiber;
}


G4Material* PHG4BarrelEcalDetector::GetSciGlass()
{
  static string matname = "sciglass";
  G4Material* sciglass = G4Material::GetMaterial(matname, false);  // false suppresses warning that material does not exist
  if (!sciglass)
  {
    G4double density;
    G4int ncomponents;
    sciglass = new G4Material(matname, density = 4.22 * g / cm3, ncomponents = 4, kStateSolid);
    sciglass->AddElement(G4Element::GetElement("Ba"), 0.3875);
    sciglass->AddElement(G4Element::GetElement("Gd"), 0.2146);
    sciglass->AddElement(G4Element::GetElement("Si"), 0.1369);
    sciglass->AddElement(G4Element::GetElement("O"), 0.2610);
  }
  return sciglass;
}

G4GenericTrap* PHG4BarrelEcalDetector::GetTowerTrap(std::map<std::string, towerposition>::iterator iterator)
{
  G4double zheight = iterator->second.size_height;
  std::vector<G4TwoVector> trapcoords(8);
  trapcoords[0] = G4TwoVector(-iterator->second.size_xin, -iterator->second.size_xin);
  trapcoords[1] = G4TwoVector(-iterator->second.size_xin, iterator->second.size_xin);
  trapcoords[2] = G4TwoVector(iterator->second.size_xin, iterator->second.size_xin);
  trapcoords[3] = G4TwoVector(iterator->second.size_xin, -iterator->second.size_xin);

  trapcoords[4] = G4TwoVector(-iterator->second.size_xout, -iterator->second.size_xout);
  trapcoords[5] = G4TwoVector(-iterator->second.size_xout, iterator->second.size_xout);
  trapcoords[6] = G4TwoVector(iterator->second.size_xout, iterator->second.size_xout);
  trapcoords[7] = G4TwoVector(iterator->second.size_xout, -iterator->second.size_xout);

  G4GenericTrap* block_tower = new G4GenericTrap("solid_tower",
    zheight / 2, trapcoords
  );

  return block_tower;
}

G4GenericTrap* PHG4BarrelEcalDetector::GetGlassTrap(std::map<std::string, towerposition>::iterator iterator, bool forSubtraction = false)
{
  G4double carbon_wall = m_Params->get_double_param("thickness_wall") * cm;
  G4double zheight = iterator->second.size_height;
  if(forSubtraction) zheight+=0.01*cm;
  std::vector<G4TwoVector> trapcoords(8);
  trapcoords[0] = G4TwoVector(-(iterator->second.size_xin-carbon_wall), -(iterator->second.size_xin-carbon_wall));
  trapcoords[1] = G4TwoVector(-(iterator->second.size_xin-carbon_wall), (iterator->second.size_xin-carbon_wall));
  trapcoords[2] = G4TwoVector((iterator->second.size_xin-carbon_wall), (iterator->second.size_xin-carbon_wall));
  trapcoords[3] = G4TwoVector((iterator->second.size_xin-carbon_wall), -(iterator->second.size_xin-carbon_wall));

  trapcoords[4] = G4TwoVector(-(iterator->second.size_xout-carbon_wall), -(iterator->second.size_xout-carbon_wall));
  trapcoords[5] = G4TwoVector(-(iterator->second.size_xout-carbon_wall), (iterator->second.size_xout-carbon_wall));
  trapcoords[6] = G4TwoVector((iterator->second.size_xout-carbon_wall), (iterator->second.size_xout-carbon_wall));
  trapcoords[7] = G4TwoVector((iterator->second.size_xout-carbon_wall), -(iterator->second.size_xout-carbon_wall));

  G4GenericTrap* block_tower = new G4GenericTrap("solid_tower_glass",
    zheight / 2, trapcoords
  );

  return block_tower;
}


int PHG4BarrelEcalDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  std::ifstream istream_mapping;
  istream_mapping.open(m_Params->get_string_param("mapping_file"));
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in PHG4BarrelEcalDetector: Failed to open mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
    gSystem->Exit(1);
  }

  /* loop over lines in file */
  std::string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {
    std::istringstream iss(line_mapping);

    if (line_mapping.find("BECALtower ") != string::npos)
    {
      unsigned idphi_j, ideta_k;
      G4double cx, cy, cz;
      G4double rot_z, rot_y, rot_x;
      G4double size_height, size_xin, size_xout;
      std::string dummys;
      // cout << "BECALtower " << itow << " " << 0 << " " << Lin << " " << Lout << " " << height << " " << xgrav << " " << ygrav << " " << zgrav << " " << theta0+theta1 << endl;

      if (!(iss >> dummys >> ideta_k >> idphi_j >> size_xin >> size_xout >> size_height >> cx >> cy >> cz >> rot_x >> rot_y >> rot_z))
      {
        std::cout << "ERROR in PHG4BarrelEcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      std::ostringstream towername;
      towername.str("");
      towername << m_TowerLogicNamePrefix << "_j_" << idphi_j << "_k_" << ideta_k;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.size_xin = size_xin * cm;
      tower_new.size_xout = size_xout * cm;
      tower_new.size_height = size_height * cm;
      tower_new.sizey2 = 0 * cm;
      tower_new.sizez = 0 * cm;
      tower_new.pTheta = 0;
      tower_new.centerx = cx * cm;
      tower_new.centery = cy * cm;
      tower_new.centerz = cz * cm;
      tower_new.rotx = rot_x;
      tower_new.roty = rot_y;
      tower_new.rotz = rot_z;
      tower_new.idx_j = idphi_j;
      tower_new.idx_k = ideta_k;
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
        cout << "ERROR in PHG4CrystalCalorimeterDetector: Failed to read line in mapping file " << m_Params->get_string_param("mappingtower") << endl;
        gSystem->Exit(1);
      }

      m_GlobalParameterMap.insert(make_pair(parname, parval));

      /* Update member variables for global parameters based on parsed parameter file */

      std::map<string, G4double>::iterator parit;

      parit = m_GlobalParameterMap.find("CenterZ_Shift");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("CenterZ_Shift", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("radius");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("radius", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("tower_length");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("tower_length", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("becal_length");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("becal_length", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("cone1_h");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("cone1_h", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("cone1_dz");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("cone1_dz", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("cone2_h");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("cone2_dz", parit->second);  // in cm
      }
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("cone2_dz", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("thickness_wall");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("thickness_wall", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("silicon_width_half");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("silicon_width_half", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("kapton_width_half");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("kapton_width_half", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("SIO2_width_half");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("SIO2_width_half", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("Carbon_width_half");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("Carbon_width_half", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("support_length");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("support_length", parit->second);  // in cm
      }
      parit = m_GlobalParameterMap.find("useLeadGlass");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_int_param("useLeadGlass", parit->second);
        m_useLeadGlass = m_Params->get_int_param("useLeadGlass");
        if(m_useLeadGlass) std::cout << "using lead glass as tower material" << std::endl;
      }
    }
  }

  return 0;
}
