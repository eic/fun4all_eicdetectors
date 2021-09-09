#ifndef EICCALORECO_RAWCLUSTERBUILDERKV3_H
#define EICCALORECO_RAWCLUSTERBUILDERKV3_H

#include <calobase/RawTower.h>

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;

// Basic struct for towers in clusterizer
typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_iL;
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

  RawClusterContainer *_clusters;

  float _seed_e;
  float _agg_e;
  int chkenergyconservation;
  float aggregation_margin_V3 = 0.03;

  std::string detector;
  std::string ClusterNodeName;
};

#endif
