#include "EICG4ZDCRawTowerBuilderByHitIndex.h"
#include "EICG4ZDCdetid.h"

#include <eiczdcbase/RawTowerZDCContainer.h>
#include <eiczdcbase/RawTowerZDCv1.h>

#include <eiczdcbase/RawTowerZDCGeom.h>             // for RawTowerGeom
#include <eiczdcbase/RawTowerZDCGeomv1.h>
#include <eiczdcbase/RawTowerZDCGeomContainer.h>    // for RawTowerGeomContainer

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                      // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                    // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

#include <TRotation.h>
#include <TVector3.h>

#include <cstdlib>                            // for exit
#include <exception>                           // for exception
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>                             // for pair, make_pair
#include <bitset>

using namespace std;

EICG4ZDCRawTowerBuilderByHitIndex::EICG4ZDCRawTowerBuilderByHitIndex(const std::string &name)
  : SubsysReco(name)
  , m_Towers(nullptr)
  , m_Geoms(nullptr)
  , m_Detector("NONE")
  , m_SubDetector("NONE")
  , m_MappingTowerFile("default.txt")
  , m_CaloId(RawTowerZDCDefs::NONE)
  , m_Emin(1e-9)
  , m_TowerDepth(0)
  , m_ThicknessAbsorber(0)
  , m_ThicknessScintilator(0)
{
}

int EICG4ZDCRawTowerBuilderByHitIndex::InitRun(PHCompositeNode *topNode)
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

int EICG4ZDCRawTowerBuilderByHitIndex::process_event(PHCompositeNode *topNode)
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
    if (g4hit_i->get_hit_type() != m_SubDetID) continue;

    /* encode CaloTowerID from i (xid), j (yid) index of tower / hit and calorimeter ID */
    int layerid = g4hit_i->get_layer();
    int towerid = -1;
    if(m_SubDetID==ZDCID::Scintillator){
      for(auto it = m_TowerIDtoLayerIDMap.begin(); it!=m_TowerIDtoLayerIDMap.end();it++)
	if(layerid >= it->second) towerid = it->first;
    }else {
      for(auto it = m_TowerIDtoLayerIDMap.begin(); it!=m_TowerIDtoLayerIDMap.end(); it++){
	if(layerid == it->second){
	  towerid = it->first;
	  break;
	}
      }
    }

    RawTowerZDCDefs::keytype calotowerid = RawTowerZDCDefs::encode_towerid_zdc(m_CaloId,
								     g4hit_i->get_index_i(),
								     g4hit_i->get_index_j(),
								     towerid
								     );
    
  
    /* add the energy to the corresponding tower */
    RawTowerZDCv1 *tower = dynamic_cast<RawTowerZDCv1 *>(m_Towers->getTower(calotowerid));
    if (!tower)
    {
      tower = new RawTowerZDCv1(calotowerid);
      tower->set_energy(0);
      m_Towers->AddTower(tower->get_id(), tower);
      if (Verbosity() > 2) 
      {
	cout<<m_SubDetector<<" : L = "<<layerid<<" , T = "<<towerid<<endl;
	    
        std::cout << "in:   " <<  g4hit_i->get_index_i() << "\t" << g4hit_i->get_index_j() << "\t" << towerid << std::endl;
        std::cout << "decoded: " <<  tower->get_bineta() << "\t" << tower->get_binphi()  << "\t" << tower->get_binl()  << std::endl;
      }
    }

    double hit_energy = 0;
    if(m_SubDetID==ZDCID::Crystal || m_SubDetID==ZDCID::Scintillator) hit_energy = g4hit_i->get_light_yield();
    else hit_energy = g4hit_i->get_edep();

    tower->add_ecell(layerid, hit_energy);
    tower->set_energy(tower->get_energy() + hit_energy);
    tower->add_eshower(g4hit_i->get_shower_id(), hit_energy);

  }

  float towerE = 0.;

  if (Verbosity())
  {
    towerE = m_Towers->getTotalEdep();
  }

  m_Towers->compress(m_Emin);
  if (Verbosity())
  {
    std::cout << "storing towers: "<< m_Towers->size() << std::endl;
    cout << "Energy lost by dropping towers with less than " << m_Emin
         << " energy, lost energy: " << towerE - m_Towers->getTotalEdep() << endl;
    m_Towers->identify();
    std::cout<<"Total energy: "<<m_Towers->getTotalEdep()<<std::endl;
    RawTowerZDCContainer::ConstRange begin_end = m_Towers->getTowers();
    RawTowerZDCContainer::ConstIterator iter;
    for (iter = begin_end.first; iter != begin_end.second; ++iter)
    {
      iter->second->identify();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int EICG4ZDCRawTowerBuilderByHitIndex::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void EICG4ZDCRawTowerBuilderByHitIndex::Detector(const std::string &d)
{
  m_Detector = d;

}

void EICG4ZDCRawTowerBuilderByHitIndex::SubDetector(const std::string &d)
{
  
  m_SubDetector = d;

  if(d == "ZDC_Crystal")
    m_SubDetID = ZDCID::Crystal;
  else if(d == "ZDC_SiPixel")
    m_SubDetID = ZDCID::SI_PIXEL;
  else if(d == "ZDC_SiPad")
    m_SubDetID = ZDCID::SI_PAD;
  else if(d == "ZDC_Sci")
    m_SubDetID = ZDCID::Scintillator;
  else
  {
    std::cout<<"Invalid ZDC subdetector name in EICG4ZDCRawTowerBuilderByHitIndex"
	     <<std::endl;
    exit(1);
  }
  
  m_CaloId = RawTowerZDCDefs::convert_name_to_caloid(d);
  
}

void EICG4ZDCRawTowerBuilderByHitIndex::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in EICG4ZDCRawTowerBuilderByHitIndex::CreateNodes");
  }

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in EICG4ZDCRawTowerBuilderByHitIndex::CreateNodes");
  }

  // Create the tower geometry node on the tree
  m_Geoms = new RawTowerZDCGeomContainer(m_CaloId);
  string NodeNameTowerGeometries = "TOWERGEOM_" + m_SubDetector;

  PHIODataNode<PHObject> *geomNode = new PHIODataNode<PHObject>(m_Geoms, NodeNameTowerGeometries, "PHObject");
  runNode->addNode(geomNode);

  // Find detector node (or create new one if not found)
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst(
      "PHCompositeNode", m_SubDetector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_SubDetector);
    dstNode->addNode(DetNode);
  }

  // Create the tower nodes on the tree
  m_Towers = new RawTowerZDCContainer(m_CaloId);
  string NodeNameTowers;
  if (m_SimTowerNodePrefix.empty())
  {
    // no prefix, consistent with older convension
    NodeNameTowers = "TOWER_" + m_SubDetector;
  }
  else
  {
    NodeNameTowers = "TOWER_" + m_SimTowerNodePrefix + "_" + m_SubDetector;
  }

  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_Towers, NodeNameTowers, "PHObject");
  DetNode->addNode(towerNode);

  return;
}

bool EICG4ZDCRawTowerBuilderByHitIndex::ReadGeometryFromTable()
{
  /* Stream to read table from file */
  ifstream istream_mapping;

  if(m_MappingTowerFile.find(m_SubDetector)==std::string::npos){
    cerr<<"ZDCTowerGeomManager::ReadGeometryFileFromTable - ERROR unmatched mapping file: SubDetector "<<m_SubDetector<<" : mapping file : "<<m_MappingTowerFile<<endl;
  }

  /* Open the datafile, if it won't open return an error */
  if (!istream_mapping.is_open())
  {
    istream_mapping.open(m_MappingTowerFile);
    if (!istream_mapping)
    {
      cerr << "ZDCTowerGeomManager::ReadGeometryFromTable - ERROR Failed to open mapping file " << m_MappingTowerFile << endl;
      exit(1);
    }
  }
  
  string line_mapping;

  double RotAroundX=0.;
  double RotAroundY=0.;
  double RotAroundZ=0.;
  double GlobalPlaceInX=0.;
  double GlobalPlaceInY=0.;
  double GlobalPlaceInZ=0.;

  m_TowerIDtoLayerIDMap.clear();

  while (getline(istream_mapping, line_mapping))
  {
    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != string::npos)
    {
      if (Verbosity() > 0)
      {
        cout << "EICG4ZDCRawTowerBuilderByHitIndex: SKIPPING line in mapping file: " << line_mapping << endl;
      }
      continue;
    }
    istringstream iss(line_mapping);

    /* If line starts with keyword Tower, add to tower positions */
    if (line_mapping.find("Tower ") != string::npos)
    {
      unsigned idx_T, idx_x, idx_y;
      double pos_x, pos_y, pos_z;
      double size_x, size_y, size_z;
      string dummys;

      /* read string- break if error */
      if (!(iss >> dummys >> idx_T >> idx_x >> idx_y >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z))
      {
        cerr << "ERROR in EICG4ZDCRawTowerBuilderByHitIndex: Failed to read line in mapping file " << m_MappingTowerFile << endl;
        exit(1);
      }

        /* Construct unique Tower ID */
      unsigned int temp_id = RawTowerZDCDefs::encode_towerid_zdc(m_CaloId, idx_x, idx_y, idx_T);

      /* Create tower geometry object */
      RawTowerZDCGeomv1 *temp_geo = new RawTowerZDCGeomv1(temp_id);
      temp_geo->set_center_x(pos_x);
      temp_geo->set_center_y(pos_y);
      temp_geo->set_center_z(pos_z);
      temp_geo->set_size_x(size_x);
      temp_geo->set_size_y(size_y);
      temp_geo->set_size_z(size_z);
  
      /* Insert this tower into position map */
      m_Geoms->add_tower_geometry(temp_geo);
    
    }

    /* If line does not start with keyword Tower, read as parameter */
    else if (line_mapping.find("iT_iL ")!=string::npos)
    {
      int towerid, layerid;
      string dummys;
      if(!(iss >> dummys >> towerid >> layerid))
      {
	cerr<<"ERROR in EICG4ZDCRawTowerBuilderByHitIndex: Failed to read line in mapping file "<< m_MappingTowerFile <<endl;
	exit(1);
      }
      m_TowerIDtoLayerIDMap.insert(make_pair(towerid,layerid));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      string parname;
      double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        cerr << "ERROR in EICG4ZDCRawTowerBuilderByHitIndex: Failed to read line in mapping file " << m_MappingTowerFile << endl;
        exit(1);
      }

      m_GlobalParameterMap.insert(make_pair(parname, parval));
      
      /* Update member variables for global parameters based on parsed parameter file */
      std::map<string, double>::iterator parit;

      parit = m_GlobalParameterMap.find("Gx0");
      if (parit != m_GlobalParameterMap.end())
        GlobalPlaceInX = parit->second;

      parit = m_GlobalParameterMap.find("Gy0");
      if (parit != m_GlobalParameterMap.end())
        GlobalPlaceInY = parit->second;

      parit = m_GlobalParameterMap.find("Gz0");
      if (parit != m_GlobalParameterMap.end())
        GlobalPlaceInZ = parit->second;

      parit = m_GlobalParameterMap.find("Grot_x");
      if (parit != m_GlobalParameterMap.end())
        RotAroundX = parit->second;

      parit = m_GlobalParameterMap.find("Grot_y");
      if (parit != m_GlobalParameterMap.end())
        RotAroundY = parit->second;

      parit = m_GlobalParameterMap.find("Grot_z");
      if (parit != m_GlobalParameterMap.end())
        RotAroundZ = parit->second;
      
    }
  }
    
  /* Correct tower geometries for global calorimter translation / rotation 
   * after reading parameters from file */
  RawTowerZDCGeomContainer::ConstRange all_towers = m_Geoms->get_tower_geometries();

  for (RawTowerZDCGeomContainer::ConstIterator it = all_towers.first;
       it != all_towers.second; ++it)
  {
    double x_temp = it->second->get_center_x();
    double y_temp = it->second->get_center_y();
    double z_temp = it->second->get_center_z();

    TVector3 v_temp_r(x_temp, y_temp, z_temp);

    /* Rotation */
    v_temp_r.RotateX(RotAroundX);
    v_temp_r.RotateY(RotAroundY);
    v_temp_r.RotateZ(RotAroundZ);

    /* Translation */
    double x_temp_rt = v_temp_r.X() + GlobalPlaceInX;
    double y_temp_rt = v_temp_r.Y() + GlobalPlaceInY;
    double z_temp_rt = v_temp_r.Z() + GlobalPlaceInZ;

    /* Update tower geometry object */
    it->second->set_center_x(x_temp_rt);
    it->second->set_center_y(y_temp_rt);
    it->second->set_center_z(z_temp_rt);

    if (Verbosity() > 2)
    {
      cout << "* Local tower x y z : " << x_temp << " " << y_temp << " " << z_temp << endl;
      cout << "* Globl tower x y z : " << x_temp_rt << " " << y_temp_rt << " " << z_temp_rt << endl;
    }
  }

  return true;
}
