#ifndef EICCALORECO_RAWCLUSTERBUILDERHELPER_H
#define EICCALORECO_RAWCLUSTERBUILDERHELPER_H

#include <calobase/RawTower.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>

#include <string>


class RawClusterBuilderHelper : public SubsysReco
{
    public:
        typedef struct {
        float tower_E;
        int tower_iEta;
        int tower_iPhi;
        int tower_trueID;
        RawTower *twr;
        } towersStrct;
        static bool towerECompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }

        RawClusterBuilderHelper(const std::string &name);
        ~RawClusterBuilderHelper() override {};

        int InitRun(PHCompositeNode *topNode) override;
----->  int process_event(PHCompositeNode *topNode) override;
        int End(PHCompositeNode *topNode) override;

        void Detector(const std::string &d) { detector = d; }
        void set_seed_e(const float e) { _seed_e = e; }
        void set_agg_e(const float e) { _agg_e = e; }
    
    protected:
        float _seed_e;
        float _agg_e;
        std::string detector;
        std::string ClusterNodeName;

        RawClusterContainer *_clusters;

----->  void cluster(std::vector<towersStrct> &input_towers, uint caloId);
        
        int caloTowersPhi(int caloID);
        bool IsForwardCalorimeter(int caloID);
        void CreateNodes(PHCompositeNode *topNode);
    

};

#endif // EICCALORECO_RAWCLUSTERBUILDERHELPER_H