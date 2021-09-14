#include <RawClusterBuilderkMA.h>
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

#include <TMath.h>

#include <algorithm>
#include <cassert>
#include <string>
#include <utility>


RawClusterBuilderkMA::RawClusterBuilderkMA(const std::string &name)
  : RawClusterBuilderHelper(name)
{
}

void RawClusterBuilderkMA::cluster(std::vector<towersStrct> &input_towers, uint caloId)
{
  std::sort(input_towers.begin(), input_towers.end(), &towerECompare);
  std::vector<int> clslabels;
  std::vector<towersStrct> cluster_towers;
  while (!input_towers.empty()) {
    cluster_towers.clear();
    clslabels.clear();
    
    // always start with highest energetic tower
    if(input_towers.at(0).tower_E > _seed_e){
      if (Verbosity() > 2) std::cout << "new cluster" << std::endl;
      // fill seed cell information into current cluster
      cluster_towers.push_back(input_towers.at(0));
      RawCluster *cluster = new RawClusterv1();
      _clusters->AddCluster(cluster);
      cluster->addTower(input_towers.at(0).twr->get_id(),input_towers.at(0).tower_E);
      clslabels.push_back(input_towers.at(0).tower_trueID);
      std::cout << "running MA" << std::endl;
      std::vector<towersStrct> cluster_towers;
      std::vector<int> clslabels;
      // remove seed tower from sample
      input_towers.erase(input_towers.begin());
      for (int tit = 0; tit < (int)cluster_towers.size(); tit++){
        // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
        int iEtaTwr = cluster_towers.at(tit).tower_iEta;
        int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
        int iLTwr   = cluster_towers.at(tit).tower_iL;
        int refC = 0;
        for (int ait = 0; ait < (int)input_towers.size(); ait++){
          int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
          int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;
          int iLTwrAgg   = input_towers.at(ait).tower_iL;
          
          if (!IsForwardCalorimeter(caloId)) {
            if (iPhiTwr < 5 && iPhiTwrAgg > caloTowersPhi(caloId)-5){
              iPhiTwrAgg= iPhiTwrAgg-caloTowersPhi(caloId);
            }
            if (iPhiTwr > caloTowersPhi(caloId)-5 && iPhiTwrAgg < 5){
              iPhiTwr= iPhiTwr-caloTowersPhi(caloId);
            }
          }
          
          int deltaL    = TMath::Abs(iLTwrAgg-iLTwr) ;
          int deltaPhi  = TMath::Abs(iPhiTwrAgg-iPhiTwr) ;
          int deltaEta  = TMath::Abs(iEtaTwrAgg-iEtaTwr) ;
          bool neighbor = (deltaL+deltaPhi+deltaEta == 1);
          bool corner2D = (deltaL == 0 && deltaPhi == 1 && deltaEta == 1) || (deltaL == 1 && deltaPhi == 0 && deltaEta == 1) || (deltaL == 1 && deltaPhi == 1 && deltaEta == 0);          
          // first condition asks for V3-like neighbors, while second condition also checks diagonally attached towers
          if(neighbor || corner2D ){

            // only aggregate towers with lower energy than current tower
            if(caloId != RawTowerDefs::LFHCAL){
              if(input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + _agg_e)) continue;
            } 
            cluster_towers.push_back(input_towers.at(ait));
            cluster->addTower(input_towers.at(ait).twr->get_id(), input_towers.at(ait).tower_E);  // Add tower to cluster)
            if(!(std::find(clslabels.begin(), clslabels.end(), input_towers.at(ait).tower_trueID) != clslabels.end())){
              clslabels.push_back(input_towers.at(ait).tower_trueID);
            }
            if (Verbosity() > 2) std::cout << "aggregated: "<< iEtaTwrAgg << "\t" << iPhiTwrAgg << "\t" << iLTwrAgg << "\t E:" << input_towers.at(ait).tower_E << "\t reference: "<< refC << "\t"<< iEtaTwr << "\t" << iPhiTwr << "\t" << iLTwr << "\t cond.: \t"<< neighbor << "\t" << corner2D << "\t  diffs: " << deltaEta << "\t" << deltaPhi << "\t" << deltaL<< std::endl;

            input_towers.erase(input_towers.begin()+ait);
            ait--;
            refC++;
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

