#include <RawClusterBuilderHelper.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterDefs.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

RawClusterBuilderHelper::RawClusterBuilderHelper(const std::string &name)
  : SubsysReco(name)
  , _seed_e(0.5)
  , _agg_e(0.1)
  , detector("NONE")
  , _clusters(nullptr)
{
}

int RawClusterBuilderHelper::InitRun(PHCompositeNode *topNode)
{
  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << PHWHERE << ": " << e.what() << std::endl;
    throw;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int RawClusterBuilderHelper::process_event(PHCompositeNode *topNode)
{
  std::string towernodename = "TOWER_CALIB_" + detector;
  // Grab the towers
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, towernodename);
  if (!towers)
  {
    std::cout << PHWHERE << ": Could not find node " << towernodename << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  std::string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << towergeomnodename << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // make the list of towers above threshold
  std::vector<towersStrct> input_towers;
  // towers in the current cluster
  int towers_added = 0;
  RawTowerContainer::ConstRange begin_end = towers->getTowers();
  for (RawTowerContainer::ConstIterator itr = begin_end.first; itr != begin_end.second; ++itr)
  {
    RawTower *tower = itr->second;
    RawTowerDefs::keytype towerid = itr->first;
    if (tower->get_energy() > _agg_e)
    {
      towersStrct tempTower;
      tempTower.tower_E = tower->get_energy();
      tempTower.tower_iEta = tower->get_bineta();
      tempTower.tower_iPhi = tower->get_binphi();
      tempTower.tower_iL = 0;
      if (towerid == RawTowerDefs::LFHCAL)
      {
        tempTower.tower_iL = tower->get_binl();
      }
      tempTower.tower_trueID = towerid;  // currently unsigned -> signed, will this matter?
      tempTower.twr = itr->second;
      input_towers.push_back(tempTower);
      towers_added++;
    }
  }

  cluster(input_towers, towers->getCalorimeterID());

  // Sum x, y, z, e
  // from https://github.com/ECCE-EIC/coresoftware/blob/ae0526adf82f49cb8906d447411b90287de6a56e/offline/packages/CaloReco/RawClusterBuilderGraph.cc#L202
  for (const auto &cluster_pair : _clusters->getClustersMap())
  {
    RawClusterDefs::keytype clusterid = cluster_pair.first;
    RawCluster *cluster = cluster_pair.second;

    assert(cluster);
    assert(cluster->get_id() == clusterid);

    double sum_x(0);
    double sum_y(0);
    double sum_z(0);
    double sum_e(0);

    for (const auto tower_pair : cluster->get_towermap())
    {
      const RawTower *rawtower = towers->getTower(tower_pair.first);
      const RawTowerGeom *rawtowergeom = towergeom->get_tower_geometry(tower_pair.first);

      assert(rawtower);
      assert(rawtowergeom);
      const double e = rawtower->get_energy();

      sum_e += e;

      if (e > 0)
      {
        sum_x += e * rawtowergeom->get_center_x();
        sum_y += e * rawtowergeom->get_center_y();
        sum_z += e * rawtowergeom->get_center_z();
      }
    }  //     for (const auto tower_pair : cluster->get_towermap())

    cluster->set_energy(sum_e);

    if (sum_e > 0)
    {
      sum_x /= sum_e;
      sum_y /= sum_e;
      sum_z /= sum_e;

      cluster->set_r(sqrt(sum_y * sum_y + sum_x * sum_x));
      cluster->set_phi(atan2(sum_y, sum_x));

      cluster->set_z(sum_z);
    }

    if (Verbosity() > 1)
    {
      std::cout << "RawClusterBuilderGraph constucted ";
      cluster->identify();
    }  //  for (const auto & cluster_pair : _clusters->getClustersMap())

    // The output of the cluster will be in the RawClusterContainer class
  }
  if (Verbosity() > 1)
  {
    std::cout << "Found " << _clusters->getClustersMap().size() << " clusters in " << Name() << std::endl;
    for (const auto &cluster_pair : _clusters->getClustersMap())
    {
      std::cout << "\n\tnTowers: " << cluster_pair.second->getNTowers() << std::endl;
      std::cout << "\tE: " << cluster_pair.second->get_energy() << "\tPhi: " << cluster_pair.second->get_phi() << std::endl;
      std::cout << "\tX: " << cluster_pair.second->get_x();
      std::cout << "\tY: " << cluster_pair.second->get_y();
      std::cout << "\tZ: " << cluster_pair.second->get_z() << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterBuilderHelper::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

bool RawClusterBuilderHelper::IsForwardCalorimeter(int caloID)
{
  switch (caloID)
  {
  case RawTowerDefs::DRCALO:
    return true;
  case RawTowerDefs::FHCAL:
    return true;
  case RawTowerDefs::FEMC:
    return true;
  case RawTowerDefs::EHCAL:
    return true;
  case RawTowerDefs::EEMC:
    return true;
  case RawTowerDefs::HCALIN:
    return false;
  case RawTowerDefs::HCALOUT:
    return false;
  case RawTowerDefs::CEMC:
    return false;
  case RawTowerDefs::EEMC_crystal:
    return true;
  case RawTowerDefs::EEMC_glass:
    return true;
  case RawTowerDefs::LFHCAL:
    return true;
  case RawTowerDefs::BECAL:
    return false;
  default:
    std::cout << "IsForwardCalorimeter: caloID " << caloID << " not defined, returning false" << std::endl;
    return false;
  }
  return false;
}

int RawClusterBuilderHelper::caloTowersPhi(int caloID)
{
  switch (caloID)
  {
  case RawTowerDefs::CEMC:
    return 100;
  case RawTowerDefs::HCALIN:
    return 64;
  case RawTowerDefs::HCALOUT:
    return 64;
  case RawTowerDefs::EEMC_glass:
    return -1;
  case RawTowerDefs::BECAL:
    return 128;
  default:
    return 0;
  }
}

void RawClusterBuilderHelper::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Grab the cEMC node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(detector);
    dstNode->addNode(DetNode);
  }

  _clusters = new RawClusterContainer();
  ClusterNodeName = "CLUSTER_" + detector;
  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName, "PHObject");
  DetNode->addNode(clusterNode);
}
