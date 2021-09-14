#include "RawClusterBuilderkV3.h"
#include "RawClusterBuilderHelper.h"

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
#include <utility>
#include <vector>

using namespace std;

RawClusterBuilderkV3::RawClusterBuilderkV3(const std::string &name)
  : RawClusterBuilderHelper(name)
{
}

void RawClusterBuilderkV3::cluster(std::vector<towersStrct> &input_towers, uint caloId)
{
  // Next we'll sort the towers from most energetic to least
  // This is from https://github.com/FriederikeBock/AnalysisSoftwareEIC/blob/642aeb13b13271820dfee59efe93380e58456289/treeAnalysis/clusterizer.cxx#L281
  std::sort(input_towers.begin(), input_towers.end(), &towerECompare);
  std::vector<towersStrct> cluster_towers;
  // And run kV3 clustering
  while (!input_towers.empty())
  {
    cluster_towers.clear();
    // always start with highest energetic tower
    if (input_towers.at(0).tower_E > _seed_e)
    {
      RawCluster *cluster = new RawClusterv1();
      _clusters->AddCluster(cluster);  // Add cluster to cluster container
      // fill seed cell information into current cluster
      cluster->addTower(input_towers.at(0).twr->get_id(), input_towers.at(0).tower_E);
      // std::cout << "Started new cluster! " << input_towers.at(0).tower_E << std::endl;
      cluster_towers.push_back(input_towers.at(0));
      // kV3 Clustering
      input_towers.erase(input_towers.begin());
      for (int tit = 0; tit < (int) cluster_towers.size(); tit++)
      {
        // Now go recursively to the next 4 neighbours and add them to the cluster if they fulfill the conditions
        int iEtaTwr = cluster_towers.at(tit).tower_iEta;
        int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
        for (int ait = 0; ait < (int) input_towers.size(); ait++)
        {
          int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
          int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;

          if (!IsForwardCalorimeter(caloId))
          {
            if (iPhiTwr < 5 && iPhiTwrAgg > caloTowersPhi(caloId) - 5)
            {
              iPhiTwrAgg = iPhiTwrAgg - caloTowersPhi(caloId);
            }
            if (iPhiTwr > caloTowersPhi(caloId) - 5 && iPhiTwrAgg < 5)
            {
              iPhiTwr = iPhiTwr - caloTowersPhi(caloId);
            }
          }
          int deltaEta = std::abs(iEtaTwrAgg - iEtaTwr);
          int deltaPhi = std::abs(iPhiTwrAgg - iPhiTwr);

          if ((deltaEta + deltaPhi) == 1)
          {
            // only aggregate towers with lower energy than current tower
            if (input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + _agg_e)) continue;
            cluster->addTower(input_towers.at(ait).twr->get_id(), input_towers.at(ait).tower_E);  // Add tower to cluster
            // std::cout << "Added a tower to the cluster! " << input_towers.at(ait).tower_E << std::endl;
            cluster_towers.push_back(input_towers.at(ait));
            input_towers.erase(input_towers.begin() + ait);
            ait--;
          }
        }
      }
    }
    else
    {
      input_towers.clear();
    }
  }
}
