#include "PHG4CrystalCalorimeterDetector.h"
#include "PHG4CrystalCalorimeterDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Element.hh>  // for G4Element
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Types.hh>        // for G4double
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
PHG4CrystalCalorimeterDetector::PHG4CrystalCalorimeterDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_SuperDetector("NONE")
  , m_Params(parameters)
  , m_DisplayAction(dynamic_cast<PHG4CrystalCalorimeterDisplayAction*>(subsys->GetDisplayAction()))
  , _towerlogicnameprefix("CrystalCalorimeterTower")
  , m_IsActive(m_Params->get_int_param("active"))
  , m_AbsorberActive(m_Params->get_int_param("absorberactive"))
{
}

//_______________________________________________________________________
int PHG4CrystalCalorimeterDetector::IsInCrystalCalorimeter(G4VPhysicalVolume* volume) const
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
void PHG4CrystalCalorimeterDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  cout << endl << endl << "Construct Me Start Construct !!!" << endl << endl;
  
  if (Verbosity() > 0)
  {
    cout << "PHG4CrystalCalorimeterDetector: Begin Construction" << endl;
  }

  if (m_Params->get_string_param("mappingtower").empty())
  {
    cout << "ERROR in PHG4CrystalCalorimeterDetector: No tower mapping file specified. Abort detector construction." << endl;
    cout << "Please run set_string_param(\"mappingtower\", std::string filename ) first." << endl;
    exit(1);
  }

  /* Read parameters for detector construction and mapping from file */
  ParseParametersFromTable();

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  recoConsts* rc = recoConsts::instance();
  /*
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));

  G4VSolid* eemc_envelope_solid = new G4Cons("eemc_envelope_solid",
                                             m_Params->get_double_param("rMin1") * cm, m_Params->get_double_param("rMax1") * cm,
                                             m_Params->get_double_param("rMin2") * cm, m_Params->get_double_param("rMax2") * cm,
                                             m_Params->get_double_param("dz") * cm / 2.,
                                             0, 2 * M_PI);

  G4LogicalVolume* eemc_envelope_log = new G4LogicalVolume(eemc_envelope_solid, WorldMaterial, G4String("eemc_envelope"), 0, 0, 0);

  GetDisplayAction()->AddVolume(eemc_envelope_log, "Envelope");
  // Define rotation attributes for envelope cone 
  G4RotationMatrix eemc_rotm;
  eemc_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  eemc_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  eemc_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);

  // Place envelope cone in simulation 
  //  ostringstream name_envelope;
  //  name_envelope.str("");
  string name_envelope = _towerlogicnameprefix + "_envelope";

  new G4PVPlacement(G4Transform3D(eemc_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm)),
                    eemc_envelope_log,
		    name_envelope,
		    logicWorld,
		    0,
		    false,
		    OverlapCheck());

  */
    
  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  /* Place calorimeter tower within envelope */
  //  PlaceTower(eemc_envelope_log, singletower);
  PlaceTower(logicWorld, singletower);

  return;
}

//_______________________________________________________________________
G4LogicalVolume* PHG4CrystalCalorimeterDetector::ConstructTower()
{

  cout << endl << endl << "Construct tower !!!" << endl << endl;
  
  if (Verbosity() > 0)
  {
    cout << "PHG4CrystalCalorimeterDetector: Build logical volume for single tower..." << endl;
  }


  /* dimesnions of full tower */
  //  G4double carbon_thickness = 0.009 * cm;
  //  G4double airgap_crystal_carbon = 0.012 * cm;
  G4double carbon_thickness = m_Params->get_double_param("carbon_gap") * cm;
  G4double airgap_crystal_carbon = m_Params->get_double_param("air_gap") * cm;
  G4double tower_dx = m_Params->get_double_param("crystal_dx") * cm + 2 * (carbon_thickness + airgap_crystal_carbon);
  G4double tower_dy = m_Params->get_double_param("crystal_dy") * cm + 2 * (carbon_thickness + airgap_crystal_carbon);
  G4double tower_dz = m_Params->get_double_param("crystal_dz") * cm;
  G4double lead_dx = tower_dx - 2.0 * carbon_thickness;
  G4double lead_dy = tower_dy - 2.0 * carbon_thickness;

  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));

  
  /* create logical volume for single tower */
  // Building the single tower mother volume first
  // Then the crystal/sci-glass will be put in this volume
  // The shell will also placed in. The rest space leave for the air gap

  
  
  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"), tower_dx / 2.0, tower_dy / 2.0, tower_dz / 2.0);
  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid, WorldMaterial, "single_tower_logic", 0, 0, 0);
  
  /* create geometry volume for crystal inside single_tower */
  G4VSolid* solid_crystal = new G4Box(G4String("single_crystal_solid"),
                                      m_Params->get_double_param("crystal_dx") * cm / 2.0,
                                      m_Params->get_double_param("crystal_dy") * cm / 2.0,
                                      m_Params->get_double_param("crystal_dz") * cm / 2.0);
  
  /* create geometry volume for frame (carbon fiber shell) inside single_tower */
  G4VSolid* Carbon_hunk_solid = new G4Box(G4String("Carbon_hunk_solid"),
                                          tower_dx / 2.0,
                                          tower_dy / 2.0,
                                          ((tower_dz / 2.0) - 1 * mm));

  
  G4VSolid* lead_solid = new G4Box(G4String("lead_solid"), lead_dx / 2.0, lead_dy / 2.0, tower_dz / 2.0);

  G4SubtractionSolid* Carbon_shell_solid = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
                                                                  Carbon_hunk_solid,
                                                                  lead_solid,
                                                                  0,
                                                                  G4ThreeVector(0.00 * mm, 0.00 * mm, 0.00 * mm));

  
  /* create logical volumes for crystal inside single_tower */  
  G4double M_para = m_Params->get_double_param("material");
  G4Material *material_Scin, *material_shell = GetCarbonFiber();

  G4Element* ele_O = new G4Element("Oxygen", "O", 8., 16.00*g/mole);
  G4Element* ele_Si = new G4Element("Silicon", "Si", 14., 28.09*g/mole);
  G4Element* ele_B = new G4Element("Boron", "B", 5., 10.811*g/mole);
  G4Element* ele_Na = new G4Element("Sodium", "Na", 11., 22.99*g/mole);
  G4Element* ele_Mg = new G4Element("Magnesium", "Mg", 12., 24.30*g/mole);
  G4Element* ele_Pb = new G4Element("Lead", "Pb", 82., 207.2*g/mole);
  G4Element* ele_Ba = new G4Element("Barium", "Ba", 56., 137.3*g/mole);
  G4Element* ele_Gd = new G4Element("Gadolinium", "Gd", 64., 157.3*g/mole);
        
  
  if( (M_para > 0.) && (M_para < 1.) )
    {
      material_Scin = G4Material::GetMaterial("G4_PbWO4");
      cout << "Set G4_PbWO4..." << endl;
    }
  else if( (M_para > 1.) && (M_para < 2.) )
    {
      material_Scin = G4Material::GetMaterial("G4_GLASS_LEAD");
      cout << "Set G4_GLASS_LEAD..." << endl;      
    }
  else if( (M_para > 2.) && (M_para < 3.) )
    {
      material_Scin = G4Material::GetMaterial("G4_BARIUM_SULFATE");
      cout << "Set G4_BARIUM_SULFATE..." << endl;      
    }
  else if( (M_para > 3.) && (M_para < 4.) )
    {
      material_Scin = G4Material::GetMaterial("G4_CESIUM_IODIDE");
      cout << "Set G4_CESIUM_IODIDE..." << endl;      
    }
  else if( (M_para > 4.) && (M_para < 5.) )
    {      
      material_Scin = new G4Material("material_Scin", 4.5*g/cm3, 5);
      material_Scin->AddElement(ele_Si, 21.9*perCent);
      material_Scin->AddElement(ele_B, 8.8*perCent);
      material_Scin->AddElement(ele_Na, 10.4*perCent);
      material_Scin->AddElement(ele_Mg, 6.5*perCent);
      material_Scin->AddElement(ele_O, 52.4*perCent);

      cout << "Set Sciglass..." << endl;
    }
  else if( (M_para > 5.) && (M_para < 6.) )
    {
      material_Scin = new G4Material("material_Scin", 9.0*g/cm3, 5);
      material_Scin->AddElement(ele_Si, 21.9*perCent);
      material_Scin->AddElement(ele_B, 8.8*perCent);
      material_Scin->AddElement(ele_Na, 10.4*perCent);
      material_Scin->AddElement(ele_Mg, 6.5*perCent);
      material_Scin->AddElement(ele_O, 52.4*perCent);

      cout << "Set heavier Sciglass..." << endl;
    }
  else if( (M_para > 6.) && (M_para < 7.) )
    {
      material_Scin = new G4Material("material_Scin", 4.5*g/cm3, 3);
      material_Scin->AddElement(ele_Si, 21.9*perCent);
      material_Scin->AddElement(ele_O, 52.4*perCent);
      material_Scin->AddElement(ele_Pb, 25.7*perCent);

      cout << "Set Sciglass contained lead..." << endl;
    }
  else if( (M_para > 7.) && (M_para < 8.) )
    {
      material_Scin = new G4Material("material_Scin", 4.22*g/cm3, 4);
      material_Scin->AddElement(ele_O, 0.261);
      material_Scin->AddElement(ele_Ba, 0.3875);
      material_Scin->AddElement(ele_Si, 0.1369);
      material_Scin->AddElement(ele_Gd, 0.2146);

      cout << "Set Sciglass from Nathely" << endl;
    }
  else if( (M_para > 8.) && (M_para < 9.) )
    {
      material_Scin = new G4Material("material_Scin", 3.8*g/cm3, 3);
      material_Scin->AddElement(ele_O, 0.293);
      material_Scin->AddElement(ele_Ba, 0.502);
      material_Scin->AddElement(ele_Si, 0.205);

      cout << "Set Sciglass from g4e" << endl;
    }
  
  
  
  
  G4LogicalVolume* logic_crystal = new G4LogicalVolume(solid_crystal, material_Scin, "single_crystal_logic", 0, 0, 0);
  G4double Cr = m_Params->get_double_param("color_R"),
           Cg = m_Params->get_double_param("color_G"),
           Cb = m_Params->get_double_param("color_B");
  G4VisAttributes* scin_vis = new G4VisAttributes(G4Colour(Cr, Cg, Cb));
  logic_crystal->SetVisAttributes(scin_vis);
  GetDisplayAction()->AddVolume(logic_crystal, "Crystal");

  G4LogicalVolume* logic_shell = new G4LogicalVolume(Carbon_shell_solid, material_shell, "single_absorber_logic", 0, 0, 0);
  GetDisplayAction()->AddVolume(logic_shell, "CarbonShell");

  
  /* Place structural frame in logical tower volume */
  // ostringstream name_shell;
  // name_shell.str("");
  // name_shell << _towerlogicnameprefix << "_single_absorber";
  string name_shell = _towerlogicnameprefix + "_single_shell";
  G4VPhysicalVolume* physvol = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logic_shell, name_shell, single_tower_logic, 0, 0, OverlapCheck());
  m_PassiveVolumeSet.insert(physvol);

  
  /* Place crystal in logical tower volume */
  string name_crystal = _towerlogicnameprefix + "_single_crystal";
  physvol = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logic_crystal, name_crystal, single_tower_logic, 0, 0, OverlapCheck());
  m_ActiveVolumeSet.insert(physvol);
  

  
  if (Verbosity() > 0)
    cout << "PHG4CrystalCalorimeterDetector: Building logical volume for single tower done." << endl;
  
  return single_tower_logic;
}



int PHG4CrystalCalorimeterDetector::PlaceTower(G4LogicalVolume* eemcenvelope, G4LogicalVolume* singletower)
{

  //  cout << endl << endl << "Start Construct !!!" << endl << endl;
  
  /* Loop over all tower positions in vector and place tower */
  typedef std::map<std::string, towerposition>::iterator it_type;

  for (it_type iterator = _map_tower.begin(); iterator != _map_tower.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4CrystalCalorimeterDetector: Place tower " << iterator->first
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

int PHG4CrystalCalorimeterDetector::ParseParametersFromTable()
{

  //  cout << endl << endl << "Start Construct !!!" << endl << endl;
    
  /* Open the datafile, if it won't open return an error */
  ifstream istream_mapping(m_Params->get_string_param("mappingtower"));
  if (!istream_mapping.is_open())
  {
    cout << "ERROR in PHG4CrystalCalorimeterDetector: Failed to open mapping file " << m_Params->get_string_param("mappingtower") << endl;
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
        cout << "PHG4CrystalCalorimeterDetector: SKIPPING line in mapping file: " << line_mapping << endl;
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
        cout << "ERROR in PHG4CrystalCalorimeterDetector: Failed to read line in mapping file " << m_Params->get_string_param("mappingtower") << endl;
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
        cout << "ERROR in PHG4CrystalCalorimeterDetector: Failed to read line in mapping file " << m_Params->get_string_param("mappingtower") << endl;
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

  
  return 0;
}

G4Material* PHG4CrystalCalorimeterDetector::GetCarbonFiber()
{

  cout << endl << endl << "Get carbon fiber !!!" << endl << endl;
    
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
