#include "B0RawTowerBuilderByHitIndex.h"

#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerv1.h>

#include <calobase/RawTower.h>               // for RawTower
#include <calobase/RawTowerDefs.h>           // for convert_name_to_caloid
#include <calobase/RawTowerGeom.h>           // for RawTowerGeom
#include <calobase/RawTowerGeomContainer.h>  // for RawTowerGeomContainer
#include <calobase/RawTowerGeomContainerv1.h>
#include <calobase/RawTowerGeomv3.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <pdbcalbase/PdbParameterMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>

#include <TRotation.h>
#include <TSystem.h>
#include <TVector3.h>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Types.hh>        // for G4double, G4int
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>  // for pair, make_pair

B0RawTowerBuilderByHitIndex::B0RawTowerBuilderByHitIndex(const std::string &name)
  : SubsysReco(name)
  , m_Towers(nullptr)
  , m_Geoms(nullptr)
  , m_Detector("NONE")
  , m_MappingTowerFile("default.txt")
  , m_CaloId(RawTowerDefs::NONE)
  , m_GlobalPlaceInX(0)
  , m_GlobalPlaceInY(0)
  , m_GlobalPlaceInZ(0)
  , m_RotInX(0)
  , m_RotInY(0)
  , m_RotInZ(0)
  , m_Emin(1e-6)
{
}

int B0RawTowerBuilderByHitIndex::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    //exit(1);
  }

  try
  {
    ReadGeometryFromTable();
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    //exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int B0RawTowerBuilderByHitIndex::process_event(PHCompositeNode *topNode)
{
  // get hits
  std::string NodeNameHits = "G4HIT_" + m_Detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, NodeNameHits);
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << NodeNameHits << std::endl;
    exit(1);
  }

  // loop over all hits in the event
  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();

  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
  {
    PHG4Hit *g4hit_i = hiter->second;

    // Don't include hits with zero energy
    if (g4hit_i->get_edep() <= 0 && g4hit_i->get_edep() != -1) continue;
    if (g4hit_i->get_index_j() < 0 || g4hit_i->get_index_k() < 0) continue;

    /* encode CaloTowerID from j, k index of tower / hit and calorimeter ID */
    //cout << "HIT indexes :"<<g4hit_i->get_index_j()<<"\t"<<g4hit_i->get_index_k()<<"\t"<<g4hit_i->get_x(0)<<"\t"<<g4hit_i->get_y(0)<<"\t"<<g4hit_i->get_light_yield()<<endl;
    RawTowerDefs::keytype calotowerid = RawTowerDefs::encode_towerid(m_CaloId,
                                                                     g4hit_i->get_index_j(),
                                                                     g4hit_i->get_index_k());
    /* add the energy to the corresponding tower */
    RawTowerv1 *tower = dynamic_cast<RawTowerv1 *>(m_Towers->getTower(calotowerid));
    if (!tower)
    {
      tower = new RawTowerv1(calotowerid);
      tower->set_energy(0);
      m_Towers->AddTower(tower->get_id(), tower);
    }

    tower->add_ecell((g4hit_i->get_index_j() << 16) + g4hit_i->get_index_k(), g4hit_i->get_light_yield());
    tower->set_energy(tower->get_energy() + g4hit_i->get_light_yield());
    tower->add_eshower(g4hit_i->get_shower_id(), g4hit_i->get_edep());
  }

  float towerE = 0.;

  if (Verbosity())
  {
    towerE = m_Towers->getTotalEdep();
    std::cout << "towers before compression: " << m_Towers->size() << "\t" << m_Detector << std::endl;
  }
  m_Towers->compress(m_Emin);
  if (Verbosity())
  {
    std::cout << "storing towers: " << m_Towers->size() << std::endl;
    std::cout << "Energy lost by dropping towers with less than " << m_Emin
         << " energy, lost energy: " << towerE - m_Towers->getTotalEdep() << std::endl;
    m_Towers->identify();
    RawTowerContainer::ConstRange begin_end = m_Towers->getTowers();
    RawTowerContainer::ConstIterator iter;
    for (iter = begin_end.first; iter != begin_end.second; ++iter)
    {
      iter->second->identify();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int B0RawTowerBuilderByHitIndex::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void B0RawTowerBuilderByHitIndex::Detector(const std::string &d)
{
  m_Detector = d;
  m_CaloId = RawTowerDefs::convert_name_to_caloid(m_Detector);
}

void B0RawTowerBuilderByHitIndex::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in B0RawTowerBuilderByHitIndex::CreateNodes");
  }

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in B0RawTowerBuilderByHitIndex::CreateNodes");
  }

  // Create the tower geometry node on the tree
  m_Geoms = new RawTowerGeomContainerv1(RawTowerDefs::convert_name_to_caloid(m_Detector));
  std::string NodeNameTowerGeometries = "TOWERGEOM_" + m_Detector;

  PHIODataNode<PHObject> *geomNode = new PHIODataNode<PHObject>(m_Geoms, NodeNameTowerGeometries, "PHObject");
  runNode->addNode(geomNode);

  // Find detector node (or create new one if not found)
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst(
      "PHCompositeNode", m_Detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }

  // Create the tower nodes on the tree
  m_Towers = new RawTowerContainer(RawTowerDefs::convert_name_to_caloid(m_Detector));
  std::string NodeNameTowers;
  if (m_SimTowerNodePrefix.empty())
  {
    // no prefix, consistent with older convension
    NodeNameTowers = "TOWER_" + m_Detector;
  }
  else
  {
    NodeNameTowers = "TOWER_" + m_SimTowerNodePrefix + "_" + m_Detector;
  }

  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_Towers, NodeNameTowers, "PHObject");
  DetNode->addNode(towerNode);

  return;
}

bool B0RawTowerBuilderByHitIndex::ReadGeometryFromTable()
{
  /* Stream to read table from file */
  std::ifstream istream_mapping;

  /* Open the datafile, if it won't open return an error */
  if (!istream_mapping.is_open())
  {
    istream_mapping.open(m_MappingTowerFile);
    if (!istream_mapping)
    {
      std::cerr << "CaloTowerGeomManager::ReadGeometryFromTable - ERROR Failed to open mapping file " << m_MappingTowerFile << std::endl;
      exit(1);
    }
  }

  std::string line_mapping;

  while (getline(istream_mapping, line_mapping))
  {
    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != std::string::npos)
    {
      if (Verbosity() > 0)
      {
        std::cout << "B0RawTowerBuilderByHitIndex: SKIPPING line in mapping file: " << line_mapping << std::endl;
      }
      continue;
    }

    //    istringstream iss(line_mapping);

    /* If line starts with keyword Tower, add to tower positions */
    //    if (line_mapping.find("Tower ") != std::string::npos)
    //    {
    //      unsigned idx_j, idx_k, idx_l;
    //      double pos_x, pos_y, pos_z;
    //      double size_x, size_y, size_z;
    //      double rot_x, rot_y, rot_z;
    //      double type;
    //      std::string dummys;

    /* read std::string- break if error */
    //      if (!(iss >> dummys >> type >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z))
    //      {
    //        std::cerr << "ERROR in RawTowerBuilderByHitIndex: Failed to read line in mapping file " << m_MappingTowerFile << std::endl;
    //        exit(1);
    //      }

    std::istringstream iss(line_mapping);

    //    /* If line starts with keyword Tower, add to tower positions */
    //    if (line_mapping.find("Tower ") != std::string::npos)
    //    {
    unsigned idx_j, idx_k, idx_l;
    G4double pos_x, pos_y, pos_z, pos_z1;
    G4double size_x, size_y, size_z;
    //      size_x = 2 *cm;
    //      size_y = 2 * cm;
    //      size_z = 20 * cm;
    G4double rot_x, rot_y, rot_z;
    G4double dummy;
    std::string dummys;
    if (line_mapping.find("B0 ") != std::string::npos)
    {
      if (!(iss >> dummys >> pos_x >> pos_y >> pos_z >> pos_z1))
      {
        std::cout << "ERROR in B0RawTowerBuilder: Failed to read global position from  mapping file " << m_MappingTowerFile << std::endl;
        gSystem->Exit(1);
      }
      m_GlobalPlaceInX = pos_x * cm;
      m_GlobalPlaceInY = pos_y * cm;
      m_GlobalPlaceInZ = (pos_z + pos_z1) * cm;
    }
    else if (line_mapping.find("Tower ") != std::string::npos)
    {
      /* read std::string- break if error */
      if (!(iss >> dummys >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z >> dummy))
      //      if (!(iss >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z >> dummy))
      {
        std::cout << "ERROR in B0RawTowerBuilder: Failed to read line in mapping file " << m_MappingTowerFile << std::endl;
        gSystem->Exit(1);
      }

      /* Construct unique Tower ID */
      unsigned int temp_id = RawTowerDefs::encode_towerid(m_CaloId, idx_j, idx_k);
      size_x = size_x * cm;
      size_y = size_y * cm;
      size_z = size_z * cm;
      pos_x = pos_x * cm;
      pos_y = pos_y * cm;
      pos_z = 0;  //48 * cm; //!!! DO NOT move tower to B0 position wrt the B0 magnet
      /* Create tower geometry object */
      RawTowerGeom *temp_geo = new RawTowerGeomv3(temp_id);
      temp_geo->set_center_x(pos_x);
      temp_geo->set_center_y(pos_y);
      temp_geo->set_center_z(pos_z);
      temp_geo->set_size_x(size_x);
      temp_geo->set_size_y(size_y);
      temp_geo->set_size_z(size_z);
      //      temp_geo->set_tower_type((int) type);

      /* Insert this tower into position map */
      m_Geoms->add_tower_geometry(temp_geo);
    }
    else
      std::cout << "ERROR in B0RawTowerBuilder: Unknown line in Mapping File " << m_MappingTowerFile << std::endl;
    //    }
    /* If line does not start with keyword Tower, read as parameter */
    //    else
    //    {
    /* If this line is not a comment and not a tower, save parameter as std::string / value. */
    //      std::string parname;
    //      double parval;

    /* read std::string- break if error */
    //      if (!(iss >> parname >> parval))
    //      {
    //        std::cerr << "ERROR in RawTowerBuilderByHitIndex: Failed to read line in mapping file " << m_MappingTowerFile << std::endl;
    //        exit(1);
    //      }

    //      m_GlobalParameterMap.insert(make_pair(parname, parval));
    //    }
  }

  /* Correct tower geometries for global calorimter translation / rotation 
  * after reading parameters from file */
  RawTowerGeomContainer::ConstRange all_towers = m_Geoms->get_tower_geometries();

  for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
       it != all_towers.second; ++it)
  {
    double x_temp = it->second->get_center_x();
    double y_temp = it->second->get_center_y();
    double z_temp = it->second->get_center_z();

    /* Rotation */
    TRotation rot;
    rot.RotateX(m_RotInX);
    rot.RotateY(m_RotInY);
    rot.RotateZ(m_RotInZ);

    TVector3 v_temp_r(x_temp, y_temp, z_temp);
    v_temp_r.Transform(rot);
    /* Translation */
    double x_temp_rt = v_temp_r.X() + m_GlobalPlaceInX;
    double y_temp_rt = v_temp_r.Y() + m_GlobalPlaceInY;
    double z_temp_rt = v_temp_r.Z() + m_GlobalPlaceInZ;

    /* Update tower geometry object */
    it->second->set_center_x(x_temp_rt);
    it->second->set_center_y(y_temp_rt);
    it->second->set_center_z(z_temp_rt);

    if (Verbosity() > 2)
    {
      std::cout << "* Local tower x y z : " << x_temp << " " << y_temp << " " << z_temp << std::endl;
      std::cout << "* Tower size x y z : " << it->second->get_size_x() << " " << it->second->get_size_y() << " " << it->second->get_size_z() << std::endl;
      std::cout << "* Globl tower x y z : " << x_temp_rt << " " << y_temp_rt << " " << z_temp_rt << std::endl;
    }
  }

  if (Verbosity())
  {
    std::cout << "size tower geom container:" << m_Geoms->size() << "\t" << m_Detector << std::endl;
  }
  return true;
}
