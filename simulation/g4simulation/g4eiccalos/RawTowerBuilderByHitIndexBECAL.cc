#include "RawTowerBuilderByHitIndexBECAL.h"

#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerv1.h>

#include <calobase/RawTower.h>               // for RawTower
#include <calobase/RawTowerDefs.h>           // for convert_name_to_caloid
#include <calobase/RawTowerGeom.h>           // for RawTowerGeom
#include <calobase/RawTowerGeomContainer.h>  // for RawTowerGeomContainer
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv4.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TRotation.h>
#include <TVector3.h>
#include <Geant4/G4Types.hh>            // for G4double, G4int

#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>  // for pair, make_pair

using namespace std;

RawTowerBuilderByHitIndexBECAL::RawTowerBuilderByHitIndexBECAL(const std::string &name)
  : SubsysReco(name)
  , m_Towers(nullptr)
  , m_Geoms(nullptr)
  , m_Detector("NONE")
  , m_MappingTowerFile("default.txt")
  , m_CaloId(RawTowerDefs::NONE)
  , m_Emin(1e-6)
{
}

int RawTowerBuilderByHitIndexBECAL::InitRun(PHCompositeNode *topNode)
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

int RawTowerBuilderByHitIndexBECAL::process_event(PHCompositeNode *topNode)
{
  // get hits
  string NodeNameHits = "G4HIT_" + m_Detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, NodeNameHits);
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node " << NodeNameHits << endl;
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

    /* encode CaloTowerID from j, k indexBECAL of tower / hit and calorimeter ID */
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
    cout << "Energy lost by dropping towers with less than " << m_Emin
         << " energy, lost energy: " << towerE - m_Towers->getTotalEdep() << endl;
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

int RawTowerBuilderByHitIndexBECAL::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawTowerBuilderByHitIndexBECAL::Detector(const std::string &d)
{
  m_Detector = d;
  m_CaloId = RawTowerDefs::convert_name_to_caloid(m_Detector);
}

void RawTowerBuilderByHitIndexBECAL::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in RawTowerBuilderByHitIndexBECAL::CreateNodes");
  }

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in RawTowerBuilderByHitIndexBECAL::CreateNodes");
  }

  // Create the tower geometry node on the tree
  m_Geoms = new RawTowerGeomContainer_Cylinderv1(RawTowerDefs::convert_name_to_caloid(m_Detector));

  string NodeNameTowerGeometries = "TOWERGEOM_" + m_Detector;

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
  string NodeNameTowers;
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

bool RawTowerBuilderByHitIndexBECAL::ReadGeometryFromTable()
{
  /* Stream to read table from file */
  ifstream istream_mapping;

  /* Open the datafile, if it won't open return an error */
  if (!istream_mapping.is_open())
  {
    istream_mapping.open(m_MappingTowerFile);
    if (!istream_mapping)
    {
      cerr << "CaloTowerGeomManager::ReadGeometryFromTable - ERROR Failed to open mapping file " << m_MappingTowerFile << endl;
      exit(1);
    }
  }

  string line_mapping;

  while (getline(istream_mapping, line_mapping))
  {
    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != string::npos)
    {
      if (Verbosity() > 0)
      {
        cout << "RawTowerBuilderByHitIndexBECAL: SKIPPING line in mapping file: " << line_mapping << endl;
      }
      continue;
    }

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
        std::cout << "ERROR in RawTowerBuilderByHitIndexBECAL: Failed to read line in mapping file " << m_MappingTowerFile << std::endl;
        exit(1);
      }

      /* Construct unique Tower ID */
      unsigned int temp_id = RawTowerDefs::encode_towerid(m_CaloId, ideta_k, idphi_j);

      /* Create tower geometry object */
      RawTowerGeom *temp_geo = new RawTowerGeomv4(temp_id);
      temp_geo->set_center_x(cx);
      temp_geo->set_center_y(cy);
      temp_geo->set_center_z(cz);
      temp_geo->set_roty(rot_x);
      temp_geo->set_roty(rot_y);
      temp_geo->set_rotz(rot_z);

      m_Geoms->add_tower_geometry(temp_geo);
    }
    else
    {
      string parname;
      double parval;
      if (!(iss >> parname >> parval))
      {
        cout << "ERROR in PHG4BarrelCalorimeterDetector: Failed to read line in mapping file " << endl;
      }

      m_GlobalParameterMap.insert(make_pair(parname, parval));

      std::map<string, double>::iterator parit;

      parit = m_GlobalParameterMap.find("thickness_wall");
      if (parit != m_GlobalParameterMap.end())
      {
        thickness_wall = parit->second;  // in cm
      }

      parit = m_GlobalParameterMap.find("radius");
      if (parit != m_GlobalParameterMap.end())
      {
        radius = parit->second;  // in cm
      }

      parit = m_GlobalParameterMap.find("tower_length");
      if (parit != m_GlobalParameterMap.end())
      {
        tower_length = parit->second;  // in cm
      }
    }
  }

  m_Geoms->set_radius(radius + 0.5 * tower_length);
  //cout << PHWHERE << "BECAL radius set to " <<  m_Geoms->get_radius() << " cm" << endl;

  // RawTowerGeomContainer::ConstRange all_towers = m_Geoms->get_tower_geometries();

  // for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
  //      it != all_towers.second; ++it)
  // {
  //   double x_temp = it->second->get_center_x();
  //   double y_temp = it->second->get_center_y();
  //   double z_temp = it->second->get_center_z();
  //   double roty = it->second->get_roty();
  //   double rotz = it->second->get_rotz();

    // double x_tempfinal = x_temp + thickness_wall / 2 * abs(cos(roty - M_PI_2)) * cos(rotz);
    // double y_tempfinal = y_temp + thickness_wall / 2 * abs(cos(roty - M_PI_2)) * sin(rotz);
    // double z_tempfinal = z_temp + thickness_wall * sin(roty - M_PI_2);

    // it->second->set_center_x(x_tempfinal);
    // it->second->set_center_y(y_tempfinal);
    // it->second->set_center_z(z_tempfinal);

    // if (Verbosity() > 2)
    // {
    //   cout << "*** Tower x y z : " << x_temp << " " << y_temp << " " << z_temp << endl;
    // }
  // }
  if (Verbosity())
  {
    cout << "size tower geom container:" << m_Geoms->size() << "\t" << m_Detector << endl;
  }

  return true;
}
