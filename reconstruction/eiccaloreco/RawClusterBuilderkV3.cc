#include "RawClusterBuilderkV3.h"



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

using namespace std;

RawClusterBuilderkV3::RawClusterBuilderkV3(const std::string &name)
  : SubsysReco(name)
  , _clusters(nullptr)
  , _seed_e(0.0)
  , _agg_e(0.0)
  , chkenergyconservation(0)
  , detector("NONE")
{
}

int RawClusterBuilderkV3::InitRun(PHCompositeNode *topNode)
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





int RawClusterBuilderkV3::process_event(PHCompositeNode *topNode)
{
  string towernodename = "TOWER_CALIB_" + detector;
  // Grab the towers
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, towernodename);
  if (!towers) {
    std::cout << PHWHERE << ": Could not find node " << towernodename << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!towergeom) {
    cout << PHWHERE << ": Could not find node " << towergeomnodename << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // make the list of towers above threshold
  std::vector<towersStrct> input_towers;  
  // towers in the current cluster
  std::vector<towersStrct> cluster_towers;
  if (towers->getCalorimeterID() == 12) { // BECAL calo_id = 12
    RawTowerContainer::ConstRange begin_end = towers->getTowers();
    for (RawTowerContainer::ConstIterator itr = begin_end.first; itr != begin_end.second; ++itr) {
      RawTower *tower = itr->second;
      RawTowerDefs::keytype towerid = itr->first;
      if (tower->get_energy() > _agg_e) {   // TODO where are aggE and E_Scaling going to be set?
        towersStrct tempTower;
        tempTower.tower_E = tower->get_energy();
        tempTower.tower_iEta = tower->get_bineta();
        tempTower.tower_iPhi = tower->get_binphi();
        tempTower.tower_trueID = towerid; // currently unsigned -> signed, will this matter?
        tempTower.twr = itr->second;
        input_towers.push_back(tempTower);
      }
    }
  }

  // Next we'll sort the towers from most energetic to least
  // This is straight from https://github.com/FriederikeBock/AnalysisSoftwareEIC/blob/642aeb13b13271820dfee59efe93380e58456289/treeAnalysis/clusterizer.cxx#L281

  RawCluster *cluster = new RawClusterv1();
  _clusters->AddCluster(cluster);

  std::sort(input_towers.begin(), input_towers.end(), &acompare);
  std::vector<int> clslabels;
  while (!input_towers.empty()) {
    cluster_towers.clear();
    clslabels.clear();
    // always start with highest energetic tower
    if(input_towers.at(0).tower_E > _seed_e){
      // fill seed cell information into current cluster
      cluster->addTower(input_towers.at(0).twr->get_id(), input_towers.at(0).tower_E);
      // kV3 Clustering
      input_towers.erase(input_towers.begin());
      for (int tit = 0; tit < (int)cluster_towers.size(); tit++){
        // Now go recursively to the next 4 neighbours and add them to the cluster if they fulfill the conditions
        int iEtaTwr = cluster_towers.at(tit).tower_iEta;
        int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
        for (int ait = 0; ait < (int)input_towers.size(); ait++){
          int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
          int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;
          
          if (iPhiTwr < 5 && iPhiTwrAgg > 128-5){  // _caloTowersPhi is the geometry? 128 for BECAL
            iPhiTwrAgg= iPhiTwrAgg-128;
          }
          if (iPhiTwr > 128-5 && iPhiTwrAgg < 5){
            iPhiTwr= iPhiTwr-128;
          }
          int deltaEta = std::abs(iEtaTwrAgg-iEtaTwr);
          int deltaPhi = std::abs(iPhiTwrAgg-iPhiTwr);

          if( (deltaEta+deltaPhi) == 1){
            // only aggregate towers with lower energy than current tower
            if(input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + aggregation_margin_V3)) continue;
            float sum_e = cluster->get_energy();
            sum_e += input_towers.at(ait).tower_E;
            cluster->set_energy(sum_e);
            cluster->addTower(input_towers.at(ait).twr->get_id(), input_towers.at(ait).tower_E);
            cluster_towers.push_back(input_towers.at(ait));
            if(!(std::find(clslabels.begin(), clslabels.end(), input_towers.at(ait).tower_trueID) != clslabels.end())){
              clslabels.push_back(input_towers.at(ait).tower_trueID);
            }
            input_towers.erase(input_towers.begin()+ait);
            ait--;
          }
        }
      }
    }

    // TODO
    // Sum x, y, z, deal with their geometry
  



  // The output of the cluster will be in the RawClusterContainer class
  }

  
  return Fun4AllReturnCodes::EVENT_OK;
}






int RawClusterBuilderkV3::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderkV3::CreateNodes(PHCompositeNode *topNode)
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
