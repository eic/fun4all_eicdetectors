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

  PlaceTower(logicWorld);
  
  return;
}

int PHG4BarrelEcalDetector::PlaceTower(G4LogicalVolume* sec)
{
  int isprojective = m_Params->get_int_param("projective");
  G4double CenterZ_Shift = 0.0;
  if(!isprojective) CenterZ_Shift = m_Params->get_double_param("CenterZ_Shift") * cm;
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

    new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(iterator->second.centerx, iterator->second.centery, iterator->second.centerz-CenterZ_Shift)),
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
  m_DisplayAction->AddVolume(shell_logic, "Carbon1");
  }
  else
  {
    m_DisplayAction->AddVolume(glass_logic, "Block2");
  m_DisplayAction->AddVolume(shell_logic, "Carbon2");
  }
  // m_DisplayAction->AddVolume(shell_logic, "Carbon");

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
    G4Element* ele_Ba = new G4Element("Barium", "Ba", 56., 137.3 * g / mole);
    G4Element* ele_Gd = new G4Element("Gadolinium", "Gd", 64., 157.3 * g / mole);
    G4double density;
    G4int ncomponents;
    sciglass = new G4Material(matname, density = 4.22 * g / cm3, ncomponents = 4, kStateSolid);
    sciglass->AddElement(ele_Ba, 0.3875);
    sciglass->AddElement(ele_Gd, 0.2146);
    // sciglass->AddElement(G4Element::GetElement("Ba"), 0.3875);
    // sciglass->AddElement(G4Element::GetElement("Gd"), 0.2146);
    sciglass->AddElement(G4Element::GetElement("Si"), 0.1369);
    sciglass->AddElement(G4Element::GetElement("O"), 0.2610);
  }
  return sciglass;
}

G4GenericTrap* PHG4BarrelEcalDetector::GetTowerTrap(std::map<std::string, towerposition>::iterator iterator)
{
  G4double zheight = iterator->second.size_height;
  std::vector<G4TwoVector> trapcoords(8);
  G4double margin = m_Params->get_double_param("margin") * cm;
  G4double x_inner = iterator->second.size_xin-margin;
  G4double x_inner_long = iterator->second.size_xinl-margin;
  G4double x_outer = iterator->second.size_xout-margin;
  G4double x_outer_long = iterator->second.size_xoutl-margin;

  int etaFlip = iterator->second.etaFlip;

  int isprojective = m_Params->get_int_param("projective");
  if(isprojective){
    if(!etaFlip){
      trapcoords[0] = G4TwoVector(-x_inner_long, -x_inner);
      trapcoords[1] = G4TwoVector(-x_inner, x_inner);
      trapcoords[2] = G4TwoVector(x_inner, x_inner);
      trapcoords[3] = G4TwoVector(x_inner_long, -x_inner);

      trapcoords[4] = G4TwoVector(-x_outer_long, -x_outer);
      trapcoords[5] = G4TwoVector(-x_outer, x_outer);
      trapcoords[6] = G4TwoVector(x_outer, x_outer);
      trapcoords[7] = G4TwoVector(x_outer_long, -x_outer);
    } else {
      trapcoords[0] = G4TwoVector(-x_inner, -x_inner);
      trapcoords[1] = G4TwoVector(-x_inner_long, x_inner);
      trapcoords[2] = G4TwoVector(x_inner_long, x_inner);
      trapcoords[3] = G4TwoVector(x_inner, -x_inner);
      trapcoords[4] = G4TwoVector(-x_outer, -x_outer);
      trapcoords[5] = G4TwoVector(-x_outer_long, x_outer);
      trapcoords[6] = G4TwoVector(x_outer_long, x_outer);
      trapcoords[7] = G4TwoVector(x_outer, -x_outer);
    }
  } else {
    trapcoords[0] = G4TwoVector(-x_inner_long, -x_inner);
    trapcoords[1] = G4TwoVector(-x_inner, x_inner);
    trapcoords[2] = G4TwoVector(x_inner, x_inner);
    trapcoords[3] = G4TwoVector(x_inner_long, -x_inner);

    trapcoords[4] = G4TwoVector(-x_outer_long, -x_outer);
    trapcoords[5] = G4TwoVector(-x_outer_long, x_outer);
    trapcoords[6] = G4TwoVector(x_outer_long, x_outer);
    trapcoords[7] = G4TwoVector(x_outer_long, -x_outer);
  }
  G4GenericTrap* block_tower = new G4GenericTrap("solid_tower",
    zheight / 2, trapcoords
  );

  return block_tower;
}

G4GenericTrap* PHG4BarrelEcalDetector::GetGlassTrap(std::map<std::string, towerposition>::iterator iterator, bool forSubtraction = false)
{
  G4double carbon_wall = m_Params->get_double_param("thickness_wall") * cm;
  G4double zheight = iterator->second.size_height;
  G4double x_inner = iterator->second.size_xin-carbon_wall;
  G4double x_outer = iterator->second.size_xout-carbon_wall;
  G4double x_inner_long = iterator->second.size_xinl-carbon_wall;
  G4double x_outer_long = iterator->second.size_xoutl-carbon_wall;
  if(forSubtraction){
    zheight+=0.1/50 * cm;
    x_inner+=0.1/50 * cm;
    x_outer+=0.1/50 * cm;
    x_inner_long+=0.1/50 * cm;
    x_outer_long+=0.1/50 * cm;
  }
  std::vector<G4TwoVector> trapcoords(8);
  int etaFlip = iterator->second.etaFlip;

  int isprojective = m_Params->get_int_param("projective");
  if(isprojective){
    if(!etaFlip){
      trapcoords[0] = G4TwoVector(-x_inner_long, -x_inner);
      trapcoords[1] = G4TwoVector(-x_inner, x_inner);
      trapcoords[2] = G4TwoVector(x_inner, x_inner);
      trapcoords[3] = G4TwoVector(x_inner_long, -x_inner);

      trapcoords[4] = G4TwoVector(-x_outer_long, -x_outer);
      trapcoords[5] = G4TwoVector(-x_outer, x_outer);
      trapcoords[6] = G4TwoVector(x_outer, x_outer);
      trapcoords[7] = G4TwoVector(x_outer_long, -x_outer);
    } else {
      trapcoords[0] = G4TwoVector(-x_inner, -x_inner);
      trapcoords[1] = G4TwoVector(-x_inner_long, x_inner);
      trapcoords[2] = G4TwoVector(x_inner_long, x_inner);
      trapcoords[3] = G4TwoVector(x_inner, -x_inner);
      trapcoords[4] = G4TwoVector(-x_outer, -x_outer);
      trapcoords[5] = G4TwoVector(-x_outer_long, x_outer);
      trapcoords[6] = G4TwoVector(x_outer_long, x_outer);
      trapcoords[7] = G4TwoVector(x_outer, -x_outer);
    }
  } else {
    trapcoords[0] = G4TwoVector(-x_inner_long, -x_inner);
    trapcoords[1] = G4TwoVector(-x_inner, x_inner);
    trapcoords[2] = G4TwoVector(x_inner, x_inner);
    trapcoords[3] = G4TwoVector(x_inner_long, -x_inner);

    trapcoords[4] = G4TwoVector(-x_outer_long, -x_outer);
    trapcoords[5] = G4TwoVector(-x_outer_long, x_outer);
    trapcoords[6] = G4TwoVector(x_outer_long, x_outer);
    trapcoords[7] = G4TwoVector(x_outer_long, -x_outer);
  }
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
      unsigned idphi_j, ideta_k, etaFlip;
      G4double cx, cy, cz;
      G4double rot_z, rot_y, rot_x;
      G4double size_height, size_xin, size_xout, size_xinl, size_xoutl;
      std::string dummys;
      // cout << "BECALtower " << itow << " " << 0 << " " << Lin << " " << Lout << " " << height << " " << xgrav << " " << ygrav << " " << zgrav << " " << theta0+theta1 << endl;

      if (!(iss >> dummys >> ideta_k >> idphi_j >> size_xin >> size_xinl >> size_xout >> size_xoutl >> size_height >> cx >> cy >> cz >> rot_x >> rot_y >> rot_z >> etaFlip))
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
      tower_new.size_xinl = size_xinl * cm;
      tower_new.size_xout = size_xout * cm;
      tower_new.size_xoutl = size_xoutl * cm;
      tower_new.size_height = size_height * cm;
      tower_new.sizey2 = 0 * cm;
      tower_new.sizez = 0 * cm;
      tower_new.etaFlip = etaFlip;
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
      parit = m_GlobalParameterMap.find("margin");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_double_param("margin", parit->second);  // in cm
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
      parit = m_GlobalParameterMap.find("projective");
      if (parit != m_GlobalParameterMap.end())
      {
        m_Params->set_int_param("projective", parit->second);
      }
    }
  }

  return 0;
}
