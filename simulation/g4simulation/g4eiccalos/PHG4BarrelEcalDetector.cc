#include "PHG4BarrelEcalDetector.h"
#include "PHG4BarrelEcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>
#include <phool/recoConsts.h>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4DisplacedSolid.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MaterialPropertiesTable.hh>  // for G4MaterialProperties...
#include <Geant4/G4MaterialPropertyVector.hh>   // for G4MaterialPropertyVector
#include <Geant4/G4SubtractionSolid.hh>
#include <TSystem.h>

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
{
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

  Radius = m_Params->get_double_param("radius")*cm;
  tower_length = m_Params->get_double_param("tower_length")*cm;

  double Length = becal_length;
  double max_radius = Radius + tower_length + elec_length + support_length;
  double pos_x1 = 0*cm;
  double pos_y1 = 0*cm;
  double pos_z1 = m_Params->get_double_param("CenterZ_Shift")*cm;

 
  G4Tubs *cylinder_solid = new G4Tubs("BCAL_SOLID",
                                       Radius, max_radius,
                                       Length/ 2.0, 0, 2*M_PI);


  G4Material *cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);


  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid, cylinder_mat,
                                       "BCAL_SOLID", 0, 0, 0);

  m_DisplayAction->AddVolume(cylinder_logic, "BCalCylinder");

  //cylinder_physi = 

  std::string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  new G4PVPlacement(0, G4ThreeVector(pos_x1, pos_y1, pos_z1), cylinder_logic, name_envelope,
                                     logicWorld, false, 0, OverlapCheck());


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
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k<< std::endl;
                //<< " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << std::endl;
    }

    int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;
  
    G4LogicalVolume*  block_logic = ConstructTower(iterator);

    m_ScintiLogicalVolSet.insert(block_logic);

    m_DisplayAction->AddVolume(block_logic, iterator->first);
    
    if(iterator->second.idx_k%2 == 0)  {
      m_DisplayAction->AddVolume(block_logic, "Block1");
    }
    else {
       m_DisplayAction->AddVolume(block_logic, "Block2");
    }

    G4LogicalVolume*  glass_logic = ConstructGlass(iterator);
    m_DisplayAction->AddVolume(glass_logic, "Glass");

    G4RotationMatrix becal_rotm;
    becal_rotm.rotateY(iterator->second.rotx);
    becal_rotm.rotateY(iterator->second.roty);
    becal_rotm.rotateZ(iterator->second.rotz);

    new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(iterator->second.centerx, iterator->second.centery, iterator->second.centerz)),
                    block_logic,
                    G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_TT")),
                    sec,
                    0, copyno, OverlapCheck());
  
    G4double posy_glass =  iterator->second.centery + th/2*abs(cos(iterator->second.roty - M_PI_2))*sin(iterator->second.rotz);
    G4double posz_glass =  iterator->second.centerz + th*sin(iterator->second.roty - M_PI_2);
    G4double posx_glass =  iterator->second.centerx + th/2*abs(cos(iterator->second.roty - M_PI_2))*cos(iterator->second.rotz);


    new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(posx_glass, posy_glass, posz_glass)),
                    glass_logic,
                    G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_Glass")),
                    sec,
                    0, copyno, OverlapCheck());

    //=================Silicon 
    G4LogicalVolume*  block_silicon = ConstructSi(iterator);
    m_DisplayAction->AddVolume(block_silicon, "Si");

    G4double pTheta   =  iterator->second.pTheta;   
    G4double theta    =  iterator->second.roty - M_PI_2;
    G4double len      = (tower_length/2) + silicon_width + overlap;
    G4double sci_sr   =  len*tan(pTheta);
    G4double sci_sz   = len;
    G4double sci_mag  = sqrt(sci_sr*sci_sr + sci_sz*sci_sz);
    

    G4double posz_si  =  iterator->second.centerz - sci_mag*sin(theta + pTheta);
    G4double posy_si  =  iterator->second.centery + sci_mag*cos(theta + pTheta)*sin(iterator->second.rotz); 
    G4double posx_si  =  iterator->second.centerx + sci_mag*cos(theta + pTheta)*cos(iterator->second.rotz);
    new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(posx_si, posy_si, posz_si)),
                    block_silicon,
                    G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_Si")),
                    sec,
                    0, copyno, OverlapCheck());

    //=================Kapton
    G4LogicalVolume*  block_kapton = ConstructKapton(iterator);
    m_DisplayAction->AddVolume(block_kapton, "Kapton");

    len               += kapton_width + overlap + silicon_width; 

    G4double kapton_sr    =  len*tan(pTheta);
    G4double kapton_sz    =  len;
    G4double kapton_mag   =  sqrt(kapton_sr*kapton_sr + kapton_sz*kapton_sz);
   
    G4double posz_kapton  =  iterator->second.centerz - kapton_mag*sin(theta + pTheta);
    G4double posy_kapton  =  iterator->second.centery + kapton_mag*cos(theta + pTheta)*sin(iterator->second.rotz); 
    G4double posx_kapton  =  iterator->second.centerx + kapton_mag*cos(theta + pTheta)*cos(iterator->second.rotz);
    new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(posx_kapton, posy_kapton, posz_kapton)),
                    block_kapton,
                    G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_Kapton")),
                    sec,
                    0, copyno, OverlapCheck());

    //=================SIO2
    G4LogicalVolume*  block_SIO2 = ConstructSIO2(iterator);
    m_DisplayAction->AddVolume(block_SIO2, "SIO2");

    len              += SIO2_width + overlap + kapton_width; 

    G4double SIO2_sr = len*tan(pTheta);
    G4double SIO2_sz = len;
    G4double SIO2_mag = sqrt(SIO2_sr*SIO2_sr + SIO2_sz*SIO2_sz);
   
    G4double posz_SIO2  =  iterator->second.centerz - SIO2_mag*sin(theta + pTheta);
    G4double posy_SIO2  =  iterator->second.centery + SIO2_mag*cos(theta + pTheta)*sin(iterator->second.rotz); 
    G4double posx_SIO2  =  iterator->second.centerx + SIO2_mag*cos(theta + pTheta)*cos(iterator->second.rotz);
    new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(posx_SIO2, posy_SIO2, posz_SIO2)),
                    block_SIO2,
                    G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("SIO2")),
                    sec,
                    0, copyno, OverlapCheck());

    //=================Carbon
    G4LogicalVolume*  block_Carbon = ConstructCarbon(iterator);
 
     if(iterator->second.idx_k%2 == 0)  {
      m_DisplayAction->AddVolume(block_Carbon, "Block1");
    }
    else {
       m_DisplayAction->AddVolume(block_Carbon, "Block2");
    }


    len              += Carbon_width + overlap + SIO2_width; 

    G4double Carbon_sr = len*tan(pTheta);
    G4double Carbon_sz = len;
    G4double Carbon_mag = sqrt(Carbon_sr*Carbon_sr + Carbon_sz*Carbon_sz);
   
    G4double posz_Carbon  =  iterator->second.centerz - Carbon_mag*sin(theta + pTheta);
    G4double posy_Carbon  =  iterator->second.centery + Carbon_mag*cos(theta + pTheta)*sin(iterator->second.rotz); 
    G4double posx_Carbon  =  iterator->second.centerx + Carbon_mag*cos(theta + pTheta)*cos(iterator->second.rotz);
    new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(posx_Carbon, posy_Carbon, posz_Carbon)),
                    block_Carbon,
                    G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("Carbon")),
                    sec,
                    0, copyno, OverlapCheck());

  }
  return 0;
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
      G4double size_z, size_y1, size_x1, size_y2, size_x2, p_Theta;
      std::string dummys;
 
      if (!(iss >> dummys >> ideta_k >> idphi_j >>  size_x1 >>  size_x2 >> size_y1 >> size_y2 >> size_z >> p_Theta >> cx >> cy >> cz >> rot_x >> rot_y >> rot_z))
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
      tower_new.sizex1  = size_x1*cm;
      tower_new.sizex2  = size_x2*cm;
      tower_new.sizey1  = size_y1*cm;
      tower_new.sizey2  = size_y2*cm;
      tower_new.sizez   = size_z*cm;
      tower_new.pTheta  = p_Theta;
      tower_new.centerx = cx*cm;
      tower_new.centery = cy*cm;
      tower_new.centerz = cz*cm;
      tower_new.roty    = rot_x;
      tower_new.roty    = rot_y;
      tower_new.rotz    = rot_z;
      tower_new.idx_j   = idphi_j;
      tower_new.idx_k   = ideta_k;
      m_TowerPostionMap.insert(make_pair(towername.str(), tower_new));

      
    } else
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

    }
  }

  return 0;
}



G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructTower(std::map<std::string, towerposition>::iterator iterator)
{

  G4Trap *block_tower = GetTowerTrap(iterator);
  G4Trap *block_glass = GetGlassTrapSubtract(iterator); 

  //*** Don't erase G4ThreeVector shift = G4ThreeVector(-th*sin(iterator->second.roty - M_PI_2), 0, th/2*abs(cos(iterator->second.roty - M_PI_2)));

  G4ThreeVector shift = G4ThreeVector(-th*sin(iterator->second.roty - M_PI_2), 0, th/2*abs(cos(iterator->second.roty - M_PI_2)));

  G4VSolid* block_solid = new G4SubtractionSolid(G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_Envelope")), block_tower, block_glass, 0, shift);

  G4Material* material_shell = GetCarbonFiber();
  assert(material_shell);
    
  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, material_shell,
                                                     G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_Tower")), 0, 0,
                                                     nullptr);
  m_ScintiLogicalVolSet.insert(block_logic);

  return block_logic;
}

G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructGlass(std::map<std::string, towerposition>::iterator iterator)
{
  G4Trap *block_solid = GetGlassTrap(iterator);   
  G4Material* material_glass = GetSciGlass();
  assert(material_glass);
  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, material_glass,
                                                     G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_Glass")), 0, 0,
                                                     nullptr);
  m_ScintiLogicalVolSet.insert(block_logic);
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
    sciglass = new G4Material(matname, density = 4.22 * g / cm3, ncomponents = 4,  kStateSolid);
    sciglass->AddElement(G4Element::GetElement("Ba"), 0.3875);
    sciglass->AddElement(G4Element::GetElement("Gd"), 0.2146);
    sciglass->AddElement(G4Element::GetElement("Si"), 0.1369);
    sciglass->AddElement(G4Element::GetElement("O"),  0.2610);
  }
  return sciglass;
}

G4Trap* PHG4BarrelEcalDetector::GetGlassTrap(std::map<std::string, towerposition>::iterator iterator)
{

  G4double size_x1 = iterator->second.sizex1/2 - th ;
  G4double size_x2 = iterator->second.sizex2/2 - th; 
  G4double size_y1 = iterator->second.sizey1/2 - th;
  G4double size_y2 = iterator->second.sizey2/2 - th; 
  G4double size_z  = iterator->second.sizez/2 - th; 

  G4Trap* block_glass = new G4Trap( "solid_glass",
      size_z,                                                                           // G4double pDz,
      iterator->second.pTheta,  0,                                                              // G4double pTheta, G4double pPhi,
      size_y1, size_x1, size_x1,      // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                        // G4double pAlp1,
      size_y2, size_x2, size_x2,      // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                         // G4double pAlp2 //
  );

  return block_glass;
}

G4Trap* PHG4BarrelEcalDetector::GetTowerTrap(std::map<std::string, towerposition>::iterator iterator)
{

  G4Trap* block_tower = new G4Trap( "solid_tower",
      iterator->second.sizez/2,                                                                           // G4double pDz,
      iterator->second.pTheta,  0,                                                              // G4double pTheta, G4double pPhi,
      iterator->second.sizey1/2, iterator->second.sizex1/2, iterator->second.sizex1/2.,      // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                        // G4double pAlp1,
      iterator->second.sizey2/2, iterator->second.sizex2/2, iterator->second.sizex2/2.,      // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                         // G4double pAlp2 //
  );

  return block_tower;
}


G4Trap* PHG4BarrelEcalDetector::GetGlassTrapSubtract(std::map<std::string, towerposition>::iterator iterator)
{

  G4double size_x1 = iterator->second.sizex1/2 - th + overlap;
  G4double size_x2 = iterator->second.sizex2/2 - th + overlap; 
  G4double size_y1 = iterator->second.sizey1/2 - th + overlap;
  G4double size_y2 = iterator->second.sizey2/2 - th + overlap; 
  G4double size_z  = iterator->second.sizez/2 - th/2  + overlap; 

  G4Trap* block_glass = new G4Trap( "solid_glass",
      size_z,                                                                           // G4double pDz,
      iterator->second.pTheta,  0,                                                              // G4double pTheta, G4double pPhi,
      size_y1, size_x1, size_x1,      // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                        // G4double pAlp1,
      size_y2, size_x2, size_x2,      // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                         // G4double pAlp2 //
  );

  return block_glass;
}

G4Trap* PHG4BarrelEcalDetector::GetSiTrap(std::map<std::string, towerposition>::iterator iterator)
{

  G4double size_x1 = iterator->second.sizex2/2 - overlap;
  G4double size_x2 = iterator->second.sizex2/2 - overlap; 
  G4double size_y1 = iterator->second.sizey2/2 - overlap;
  G4double size_y2 = iterator->second.sizey2/2 - overlap; 
  G4double size_z  = silicon_width; 

  G4Trap* block_si = new G4Trap( "Si",
      size_z,                                                                           // G4double pDz,
      iterator->second.pTheta,  0,                                                              // G4double pTheta, G4double pPhi,
      size_y1, size_x1, size_x1,      // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                        // G4double pAlp1,
      size_y2, size_x2, size_x2,      // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                         // G4double pAlp2 //
  );

  return block_si;
}

G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructSi(std::map<std::string, towerposition>::iterator iterator)
{
  G4Trap *block_solid = GetSiTrap(iterator);   
  G4Material* material_si = G4Material::GetMaterial("G4_POLYSTYRENE");
  assert(material_si);
  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, material_si,
                                                     G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_solid_Si")), 0, 0,
                                                     nullptr);
  m_ScintiLogicalVolSet.insert(block_logic);
  return block_logic;
}

G4Trap* PHG4BarrelEcalDetector::GetKaptonTrap(std::map<std::string, towerposition>::iterator iterator)
{

  G4double size_x1 = iterator->second.sizex2/2 - overlap;
  G4double size_x2 = iterator->second.sizex2/2 - overlap; 
  G4double size_y1 = iterator->second.sizey2/2 - overlap;
  G4double size_y2 = iterator->second.sizey2/2 - overlap; 
  G4double size_z  = kapton_width; 

  G4Trap* block_si = new G4Trap( "Kapton",
      size_z,                                                                           // G4double pDz,
      iterator->second.pTheta,  0,                                                              // G4double pTheta, G4double pPhi,
      size_y1, size_x1, size_x1,      // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                        // G4double pAlp1,
      size_y2, size_x2, size_x2,      // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                         // G4double pAlp2 //
  );

  return block_si;
}

G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructKapton(std::map<std::string, towerposition>::iterator iterator)
{
  G4Trap *block_solid = GetKaptonTrap(iterator);   
  G4Material* material_kapton = G4Material::GetMaterial("G4_KAPTON");
  assert(material_kapton);
  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, material_kapton,
                                                     G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_solid_Si")), 0, 0,
                                                     nullptr);
  m_ScintiLogicalVolSet.insert(block_logic);
  return block_logic;
}

G4Trap* PHG4BarrelEcalDetector::GetSIO2Trap(std::map<std::string, towerposition>::iterator iterator)
{

  G4double size_x1 = iterator->second.sizex2/2 - overlap;
  G4double size_x2 = iterator->second.sizex2/2 - overlap; 
  G4double size_y1 = iterator->second.sizey2/2 - overlap;
  G4double size_y2 = iterator->second.sizey2/2 - overlap; 
  G4double size_z  = SIO2_width; 

  G4Trap* block_SIO2 = new G4Trap( "SIO2",
      size_z,                                                                           // G4double pDz,
      iterator->second.pTheta,  0,                                                              // G4double pTheta, G4double pPhi,
      size_y1, size_x1, size_x1,      // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                        // G4double pAlp1,
      size_y2, size_x2, size_x2,      // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                         // G4double pAlp2 //
  );

  return block_SIO2;
}

G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructSIO2(std::map<std::string, towerposition>::iterator iterator)
{
  G4Trap *block_solid = GetSIO2Trap(iterator);   
  G4Material* material_SIO2 = G4Material::GetMaterial("Quartz");
  assert(material_SIO2);
  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, material_SIO2,
                                                     G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_solid_material_SIO2")), 0, 0,
                                                     nullptr);
  m_ScintiLogicalVolSet.insert(block_logic);
  return block_logic;
}


G4Trap* PHG4BarrelEcalDetector::GetCarbonTrap(std::map<std::string, towerposition>::iterator iterator)
{

  G4double size_x1 = iterator->second.sizex2/2 - overlap;
  G4double size_x2 = iterator->second.sizex2/2 - overlap; 
  G4double size_y1 = iterator->second.sizey2/2 - overlap;
  G4double size_y2 = iterator->second.sizey2/2 - overlap; 
  G4double size_z  = Carbon_width; 

  G4Trap* block_SIO2 = new G4Trap( "C",
      size_z,                                                                           // G4double pDz,
      iterator->second.pTheta,  0,                                                              // G4double pTheta, G4double pPhi,
      size_y1, size_x1, size_x1,      // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                        // G4double pAlp1,
      size_y2, size_x2, size_x2,      // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                         // G4double pAlp2 //
  );

  return block_SIO2;
}

G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructCarbon(std::map<std::string, towerposition>::iterator iterator)
{
  G4Trap *block_solid = GetCarbonTrap(iterator);   
  G4Material* material_Carbon = G4Material::GetMaterial("G4_C");
  assert(material_Carbon);
  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, material_Carbon,
                                                     G4String(string(iterator->first) + to_string(iterator->second.idx_j) +to_string(iterator->second.idx_k) + string("_solid_material_C")), 0, 0,
                                                     nullptr);
  m_ScintiLogicalVolSet.insert(block_logic);
  return block_logic;
}
