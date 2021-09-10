#include "RawDigitBuilderTTL.h"



#include <trackbase/TrkrClusterContainerv3.h>
#include <trackbase/TrkrClusterv2.h>
#include <trackbase/TrkrDefs.h>                     // for hitkey, getLayer
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssocv3.h>
// #include <calobase/RawTowerContainer.h>
// #include <calobase/RawTowerv1.h>
// #include <calobase/RawTowerv2.h>

// #include <calobase/RawTower.h>               // for RawTower
// #include <calobase/RawTowerDefs.h>           // for convert_name_to_caloid
// #include <calobase/RawTowerGeom.h>           // for RawTowerGeom
// #include <calobase/RawTowerGeomContainer.h>  // for RawTowerGeomContainer
// #include <calobase/RawTowerGeomContainerv1.h>
// #include <calobase/RawTowerGeomv3.h>


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


#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>  // for pair, make_pair

using namespace std;

RawDigitBuilderTTL::RawDigitBuilderTTL(const string &name)
  : SubsysReco(name)
  // , m_Towers(nullptr)
  // , m_Geoms(nullptr)
  // , m_hits(nullptr)
  , m_clusterlist(nullptr)
  // , m_clusterhitassoc(nullptr)
  , m_Detector("NONE")
  , m_Emin(1e-6)
{
}

int RawDigitBuilderTTL::InitRun(PHCompositeNode *topNode)
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


  std::cout << "adding TRKR_CLUSTER nodes" << std::endl;
  // Create the Cluster node if required
  auto trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    trkrclusters = new TrkrClusterContainerv3;
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
      new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawDigitBuilderTTL::process_event(PHCompositeNode *topNode)
{
  // get hits
  string NodeNameHits = "G4HIT_" + m_Detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, NodeNameHits);
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node " << NodeNameHits << endl;
    exit(1);
  }

  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // loop over all hits in the event
  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  int clusid =0;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
  {
    PHG4Hit *g4hit_i = hiter->second;

    // Don't include hits with zero energy
    if (g4hit_i->get_edep() <= 0 && g4hit_i->get_edep() != -1) continue;

    auto clus = std::make_unique<TrkrClusterv2>();

    auto ckey = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::ttl,g4hit_i->get_index_j());
    TrkrDefs::hitsetkey tmps = g4hit_i->get_strip_z_index();
    ckey |= (tmps << 8);
    tmps = g4hit_i->get_strip_y_index();
    ckey |= (tmps << 0);

    TrkrDefs::cluskey tmp = ckey;
    TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
    key |= clusid;
    clus->setClusKey(key);

    // G4double clusx,clusy,clusz;
    // GetPixelGlobalCoordinates(g4hit_i,clusx,clusy,clusz);
    clus->setPosition(0, g4hit_i->get_local_x(0));
    clus->setPosition(1, g4hit_i->get_local_y(0));
    clus->setPosition(2, g4hit_i->get_local_z(0));
    clus->setGlobal();

    m_clusterlist->addCluster(clus.release());

    clusid++;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


void RawDigitBuilderTTL::Detector(const std::string &d)
{
  m_Detector = d;
  // m_CaloId = RawTowerDefs::convert_name_to_caloid(m_Detector);
}

void RawDigitBuilderTTL::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in RawDigitBuilderTTL::CreateNodes");
  }

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in RawDigitBuilderTTL::CreateNodes");
  }

  // Find detector node (or create new one if not found)
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst(
      "PHCompositeNode", m_Detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }
  return;
}


void RawDigitBuilderTTL::PrintClusters(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1)
  {
    TrkrClusterContainer *clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (!clusterlist) return;

    cout << "================= Aftyer MvtxClusterizer::process_event() ====================" << endl;

    cout << " There are " << clusterlist->size() << " clusters recorded: " << endl;

    clusterlist->identify();

    cout << "===========================================================================" << endl;
  }

  return;
}