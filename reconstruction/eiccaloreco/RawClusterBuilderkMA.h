#ifndef EICCALORECO_RAWCLUSTERBUILDERKMA_H
#define EICCALORECO_RAWCLUSTERBUILDERKMA_H

#include <calobase/RawTower.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/SubsysReco.h>

#include <string>

typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_trueID;
  RawTower *twr;
} towersStrct2;

bool bcompare(towersStrct2 lhs, towersStrct2 rhs) { return lhs.tower_E > rhs.tower_E; }

class RawClusterBuilderkMA : public SubsysReco {
  public:
    RawClusterBuilderkMA(const std::string &name = "RawClusterBuilderkMA");
    ~RawClusterBuilderkMA() override{}

    int InitRun(PHCompositeNode *topNode) override;
    int process_event(PHCompositeNode *topNode) override;
    int End(PHCompositeNode *topNode) override;

    void Detector(const std::string &d) { detector = d; }
    void set_seed_e(const float e) { _seed_e = e; }
    void set_agg_e(const float e) { _agg_e = e; }

  private:
    void CreateNodes(PHCompositeNode *topNode);
    bool IsForwardCalorimeter(int caloID){
      switch (caloID){
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
    int caloTowersPhi(int caloID) {
      switch (caloID) {
        case RawTowerDefs::CEMC: return 100;
        case RawTowerDefs::HCALIN: return 64;
        case RawTowerDefs::HCALOUT: return 64;
        case RawTowerDefs::EEMC_glass: return -1;
        case RawTowerDefs::BECAL: return 128;
        default: return 0;
      }
    }

  RawClusterContainer *_clusters;
  float _seed_e;
  float _agg_e;
  std::string detector; 
  std::string ClusterNodeName;
};

#endif