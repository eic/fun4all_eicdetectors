#ifndef EICCALORECO_RAWCLUSTERBUILDERHELPER_H
#define EICCALORECO_RAWCLUSTERBUILDERHELPER_H

#include <calobase/RawTower.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>

#include <string>

typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_trueID;
  RawTower *twr;
} towersStrct;

// sorting function for towers
bool towerECompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }

class RawClusterBuilderHelper : public SubsysReco
{
    public:
        RawClusterBuilderHelper(const std::string &name);
        ~RawClusterBuilderHelper() override {};

        int InitRun(PHCompositeNode *topNode) override;
        int process_event(PHCompositeNode *topNode) override;
        int End(PHCompositeNode *topNode) override;

        void Detector(const std::string &d) { detector = d; }
        void set_seed_e(const float e) { _seed_e = e; }
    
    private:
        float _seed_e;
        float _agg_e;
        std::string detector;
        std::string ClusterNodeName;

        RawClusterContainer *_clusters;
        
        int caloTowersPhi(int caloID);
        bool IsForwardCalorimeter(int caloID);
        void CreateNodes(PHCompositeNode *topNode);
    

};

#endif // EICCALORECO_RAWCLUSTERBUILDERHELPER_H