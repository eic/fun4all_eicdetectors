//____________________________________________________________________________..
//
// This is a working template for the G4 Construct() method which needs to be implemented
// We wedge a method between the G4 Construct() to enable volume hierarchies on the macro
// so here it is called ConstructMe() but there is no functional difference
// Currently this installs a simple G4Box solid, creates a logical volume from it
// and places it. Put your own detector in place (just make sure all active volumes
// get inserted into the m_PhysicalVolumesSet)
//
// Rather than using hardcoded values you should consider using the parameter class
// Parameter names and defaults are set in EICG4B0Subsystem::SetDefaultParameters()
// Only parameters defined there can be used (also to override in the macro)
// to avoids typos.
// IMPORTANT: parameters have no inherent units, there is a convention (cm/deg)
// but in any case you need to multiply them here with the correct CLHEP/G4 unit
//
// The place where you put your own detector is marked with
// //begin implement your own here://
// //end implement your own here://
// Do not forget to include the G4 includes for your volumes
//____________________________________________________________________________..

#include "EICG4BwdDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
//#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Color.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <TSystem.h>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <phool/recoConsts.h> //For rc WorldMaterial

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <utility>
 
class G4VSolid;
class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
EICG4BwdDetector::EICG4BwdDetector(PHG4Subsystem *subsys,
                                 PHCompositeNode *Node,
                                 PHParameters *parameters,
                                 const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_Layer(lyr)
  , _mapping_tower_file("")
  , m_TowerLogicNamePrefix("BwdECALTower")
{
}

//_______________________________________________________________
int EICG4BwdDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  G4LogicalVolume* mylogvol = volume->GetLogicalVolume();
  if (m_LogicalVolSet.find(mylogvol) != m_LogicalVolSet.end())
  {
    return 1;
  }
  return 0;
}

int EICG4BwdDetector::GetDetId(G4VPhysicalVolume *volume) const
{
  if (IsInDetector(volume))
  {
    return 1;
  }
  return -1;
}

//_______________________________________________________________
void EICG4BwdDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  //begin implement your own here://
  // Do not forget to multiply the parameters with their respective CLHEP/G4 unit !
  if (Verbosity() > 0)
  {
    std::cout << "EICG4BwdDetector: Begin Construction" << std::endl;
  }

  cout << " !!! length = " << m_Params->get_double_param("length");
  //Print("ALL");

  G4VSolid *solidBwd = new G4Box("EICG4BwdSolid",
                                m_Params->get_double_param("width") / 2. * cm,
                                m_Params->get_double_param("height") / 2. * cm,
                                m_Params->get_double_param("length") / 2. * cm);

/*  G4VSolid *solid0 = new G4Tubs("EICG4B0ECALSolid0",
                                0.,
                                m_Params->get_double_param("outer_radius") * cm,
                                m_Params->get_double_param("length") / 2. * cm,
                                m_Params->get_double_param("startAngle") * degree,
                                m_Params->get_double_param("spanningAngle") * degree);
  G4VSolid *solidPipeHole = new G4Tubs("EICG4B0IonPipeSolid",
                                       0.,
                                       m_Params->get_double_param("pipe_hole_r") * cm,
                                       m_Params->get_double_param("length") * cm,
                                       0., 360. * degree);
  G4VSolid *solidCableHole = new G4Tubs("EICG4B0CableSolid",
                                       0.,
                                       m_Params->get_double_param("cable_hole") * cm,
                                       m_Params->get_double_param("length") * cm,
                                       0., 360. * degree);
  G4VSolid *solidPipeHole1 = new G4Box("EICG4B0PipeSolid1",
                                m_Params->get_double_param("pipe_hole") / 2. * cm,
                                m_Params->get_double_param("pipe_hole_r")  * cm,
                                m_Params->get_double_param("length")  * cm);
  G4VSolid *solid1 = new G4Tubs("EICG4B0ECALSolid1",
                                0.,
                                (m_Params->get_double_param("outer_radius") - m_Params->get_double_param("d_radius")) * cm,
                                m_Params->get_double_param("length") / 2. * cm,
                                (m_Params->get_double_param("startAngle") + m_Params->get_double_param("spanningAngle")) * degree,
                                (360 - m_Params->get_double_param("spanningAngle")) * degree);
  G4UnionSolid *solid10 = new G4UnionSolid("EICG4B0ECALSolid10", solid0, solid1);
  G4SubtractionSolid *solids = new G4SubtractionSolid("EICG4B0Solid", solid10, solidPipeHole, 0, G4ThreeVector((m_Params->get_double_param("pipe_x")+m_Params->get_double_param("pipe_hole")/2) * cm, m_Params->get_double_param("pipe_y") * cm, m_Params->get_double_param("pipe_z") * cm));
  G4SubtractionSolid *solids1 = new G4SubtractionSolid("EICG4B0Solid", solids, solidPipeHole, 0, G4ThreeVector((m_Params->get_double_param("pipe_x")-m_Params->get_double_param("pipe_hole")/2) * cm, m_Params->get_double_param("pipe_y") * cm, m_Params->get_double_param("pipe_z") * cm));
  G4SubtractionSolid *solids2 = new G4SubtractionSolid("EICG4B0Solid", solids1, solidCableHole, 0, G4ThreeVector(m_Params->get_double_param("cable_x") * cm, m_Params->get_double_param("cable_y") * cm, m_Params->get_double_param("cable_z") * cm));
  G4SubtractionSolid *solidB0 = new G4SubtractionSolid("EICG4B0Solid", solids2, solidPipeHole1, 0, G4ThreeVector(m_Params->get_double_param("pipe_x") * cm, m_Params->get_double_param("pipe_y") * cm, m_Params->get_double_param("pipe_z") * cm));*/
  G4RotationMatrix *rotm = new G4RotationMatrix();
if(m_Params->get_int_param("lightyield")){
  if (_mapping_tower_file.empty())
  {
    std::cout << "ERROR in EICG4BwdDetector: No mapping file specified. Abort detector construction." << std::endl;
//    gSystem->Exit(1);
  }
  /* Read parameters for detector construction and mappign from file */
  ParseParametersFromTable();
}
  recoConsts* rc=recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));
  G4LogicalVolume *bwd_ecal_log = new G4LogicalVolume(solidBwd,
						WorldMaterial,
                                                "BwdECAL_envelope",
						0,0,0);

  G4VisAttributes *vis = new G4VisAttributes(G4Color(0.8, 0.4, 0.2, 1.0));
  vis->SetColor(0.8, 0.4, 0.2, 1.0);
  if (m_Params->get_string_param("material") == "G4_PbWO4") vis->SetColor(0.8, 0.4, 0.2, 1.0);
  if (m_Params->get_string_param("material") == "G4_Cu") vis->SetColor(1., 0., 1., .5);
  if (m_Params->get_string_param("material") == "G4_Si") vis->SetColor(1., 1., 0., .8);
  if (m_Params->get_string_param("material") == "G4_C") vis->SetColor(1., .2, .2, .2);
  bwd_ecal_log->SetVisAttributes(vis);
  /* Place envelope cone in simulation */
  std::string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  G4VPhysicalVolume *phy = new G4PVPlacement(rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm,
                                                           m_Params->get_double_param("place_y") * cm,
                                                           m_Params->get_double_param("place_z") * cm),
                    bwd_ecal_log, name_envelope, logicWorld, 0, false, OverlapCheck());

  //Create towers for the B0 Ecal:
  /* Construct single calorimeter tower */
  if (Verbosity() > 0){
	cout << "Bwd ECal Envelope Location: x, y, z:"<<endl;
	std::cout <<"Building Calorimeter from "<<m_Params->get_string_param("material")<<endl;
  	cout<< m_Params->get_double_param("place_x")<<"\t";
	cout<< m_Params->get_double_param("place_y")<<"\t";
	cout<< m_Params->get_double_param("place_z")<<endl;
  }
if(m_Params->get_int_param("lightyield")){
  G4LogicalVolume* singletower = ConstructTower();
  /* Place calorimeter tower within envelope */
  PlaceTower(bwd_ecal_log, singletower);
}
  m_PhysicalVolumesSet.insert(phy);
  // hard code detector id to detid
  m_PhysicalVolumesDet.insert({phy, m_Params->get_double_param("detid") + 1});
  //  m_LogicalVolumesSet.insert(logical);
  //end implement your own here://

  return;
}
//_______________________________________________________________
G4LogicalVolume* EICG4BwdDetector::ConstructTower()
{
  if (Verbosity() > 0)
  {
    std::cout << "EICG4BwdDetector: Build logical volume for single tower..." << std::endl;
  }

  /* create logical volume for single tower */
  G4Material* EcalMaterial = G4Material::GetMaterial(m_Params->get_string_param("material"));
  double TowerDx = m_Params->get_double_param("tower_size") * cm;
  double TowerDy = m_Params->get_double_param("tower_size") * cm;
  double TowerDz = m_Params->get_double_param("length") * cm;
  if (Verbosity() > 0){
      std::cout << "EICG4BwdDetector: Construct tower " << m_Params->get_double_param("tower_size")<<"\t"
		<< m_Params->get_double_param("length")  << std::endl;
  }
  G4VSolid* single_tower_solid = new G4Box("single_tower_solid",
                                           TowerDx / 2.0,
                                           TowerDy / 2.0,
                                           TowerDz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            EcalMaterial,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  G4VisAttributes *vis = new G4VisAttributes(G4Color(0.8, 0.4, 0.2, .1));
  vis->SetForceSolid(false);
  single_tower_logic->SetVisAttributes(vis);
  m_LogicalVolSet.insert(single_tower_logic);
  return single_tower_logic;
}

int EICG4BwdDetector::PlaceTower(G4LogicalVolume* b0ecalenvelope, G4LogicalVolume* singletower)
{
  /* Loop over all tower positions in vector and place tower */
  for (std::map<std::string, towerposition>::iterator iterator = m_TowerPositionMap.begin(); iterator != m_TowerPositionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      std::cout << "EICG4BwdDetector: Place tower " << iterator->first
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k
                << " at x = " << iterator->second.x / cm << " , y = " << iterator->second.y / cm << " , z = " << iterator->second.z / cm << std::endl;
    }

    int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;
    new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                      singletower,
                      iterator->first,
                      b0ecalenvelope,
                      0, copyno, OverlapCheck());
//                      0, 0, OverlapCheck());
  }

  return 0;
}

int EICG4BwdDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  std::ifstream istream_mapping;
  istream_mapping.open(_mapping_tower_file);
  //istream_mapping.open("B0ECAL_mapping_v0.txt");
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in EICG4BwdDetector: Failed to open mapping file " << _mapping_tower_file << std::endl;
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
        std::cout << "EICG4BwdDetector: SKIPPING line in mapping file: " << line_mapping << std::endl;
      }
      continue;
    }

    std::istringstream iss(line_mapping);
      unsigned idx_j, idx_k, idx_l;
      G4double pos_x, pos_y, pos_z;
      double Gpos_x, Gpos_y, Gpos_z, Gpos_z1;
      G4double size_x, size_y, size_z;
      G4double rot_x, rot_y, rot_z;
      G4double dummy;
	std::string dummys;
//G4double	GlobalPlaceInX;
//G4double	GlobalPlaceInY;
//G4double	GlobalPlaceInZ;

    if (line_mapping.find("Bwd ") != std::string::npos)
    {
      if (!(iss >> dummys>> Gpos_x >> Gpos_y >> Gpos_z >> Gpos_z1 ))
      {
        std::cout << "ERROR in EICG4BwdDetector: Failed to read global position from  mapping file " << _mapping_tower_file  << std::endl;
        gSystem->Exit(1);
      }
	if (m_Params->get_double_param("global_x")!=Gpos_x || m_Params->get_double_param("global_y")!=Gpos_y ||m_Params->get_double_param("global_z")!=Gpos_z) {
	std::cout <<endl;
	std::cout << "Bwd position changed since mapping file produced or wrong mapping is used "<<_mapping_tower_file  << std::endl;
	std::cout<< m_Params->get_double_param("global_x") <<" "<< m_Params->get_double_param("global_y")<<" "<<m_Params->get_double_param("global_z") <<std::endl;
	std::cout<< Gpos_x <<" " <<Gpos_y<<" "<<Gpos_z <<std::endl;
	std::cout<< m_Params->get_double_param("global_x")-Gpos_x <<" "<< m_Params->get_double_param("global_y")-Gpos_y<<" "<<m_Params->get_double_param("global_z")-Gpos_z <<std::endl;
//	gSystem->Exit(1);
	}
//	GlobalPlaceInX=pos_x*cm;
//	GlobalPlaceInY=pos_y*cm;
//	GlobalPlaceInZ=pos_z*cm;
	}
//    /* If line starts with keyword Tower, add to tower positions */
   else if (line_mapping.find("Tower ") != std::string::npos){
      /* read string- break if error */
      if (!(iss >> dummys >>  idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z >> dummy))
      //if (!(iss >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z >> dummy))
      {
        std::cout << "ERROR in EICG4BwdDetector: Failed to read line in mapping file " << _mapping_tower_file << std::endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      std::ostringstream towername;
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
      m_TowerPositionMap.insert(make_pair(towername.str(), tower_new));
	}
	else std::cout <<"ERROR in EICG4BwdDetector: Unknown line in Mapping File " <<_mapping_tower_file << std::endl;
  }

  return 0;
}

//_______________________________________________________________
void EICG4BwdDetector::Print(const std::string &what) const
{
  std::cout << "EICG4Bwd Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}

PHParameters *EICG4BwdDetector::getParams()
{
  return m_Params;
}
