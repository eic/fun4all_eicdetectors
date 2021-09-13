#ifndef EICCALORECO_RAWCLUSTERBUILDERKV3_H
#define EICCALORECO_RAWCLUSTERBUILDERKV3_H

#include <calobase/RawTower.h>
#include <calobase/RawTowerDefs.h>

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;

// Basic struct for towers in clusterizer
typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_trueID;
  RawTower *twr;
} towersStrct;

// sorting function for towers
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }
class RawClusterBuilderkV3 : public SubsysReco
{
 public:
  RawClusterBuilderkV3(const std::string &name = "RawClusterBuilderkV3");
  ~RawClusterBuilderkV3() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void Detector(const std::string &d) { detector = d; }

  void set_seed_e(const float e) { _seed_e = e; }
  void set_agg_e(const float e) { _agg_e = e; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }

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
  int chkenergyconservation;

  std::string detector;
  std::string ClusterNodeName;
};

#endif
