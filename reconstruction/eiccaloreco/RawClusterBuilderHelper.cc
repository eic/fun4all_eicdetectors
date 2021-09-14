#include <RawClusterBuilderHelper.h>

#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawClusterDefs.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/phool.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

RawClusterBuilderHelper::RawClusterBuilderHelper(const std::string &name)
  : SubsysReco(name)
  // , _clusters(nullptr)
  // , _seed_e(0.5)
  , _agg_e(0.1)
  , detector("NONE")
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
  return Fun4AllReturnCodes::EVENT_OK;
}


int RawClusterBuilderHelper::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

bool RawClusterBuilderHelper::IsForwardCalorimeter(int caloID)
{
  switch (caloID)
  {
    case RawTowerDefs::DRCALO: return true;
    case RawTowerDefs::FHCAL: return true;
    case RawTowerDefs::FEMC: return true;
    case RawTowerDefs::EHCAL: return true;
    case RawTowerDefs::EEMC: return true;
    case RawTowerDefs::HCALIN: return false;
    case RawTowerDefs::HCALOUT: return false;
    case RawTowerDefs::CEMC: return false;
    case RawTowerDefs::EEMC_crystal: return true;
    case RawTowerDefs::EEMC_glass: return true;
    case RawTowerDefs::LFHCAL: return true;
    case RawTowerDefs::BECAL: return false;
    default:
      std::cout << "IsForwardCalorimeter: caloID " << caloID << " not defined, returning false" << std::endl;
      return false;
  }
  return false;
}

int caloTowersPhi(int caloID)
{
    switch (caloID)
    {
      case RawTowerDefs::CEMC: return 100;
      case RawTowerDefs::HCALIN: return 64;
      case RawTowerDefs::HCALOUT: return 64;
      case RawTowerDefs::EEMC_glass: return -1;
      case RawTowerDefs::BECAL: return 128;
      default: return 0;
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
