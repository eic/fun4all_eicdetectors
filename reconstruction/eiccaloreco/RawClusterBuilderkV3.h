#ifndef EICCALORECO_RAWCLUSTERBUILDERKV3_H
#define EICCALORECO_RAWCLUSTERBUILDERKV3_H

#include <RawClusterBuilderHelper.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerDefs.h>

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;

class RawClusterBuilderkV3 : public RawClusterBuilderHelper
{
  public:
    RawClusterBuilderkV3(const std::string &name = "RawClusterBuilderkV3");
    ~RawClusterBuilderkV3() override {}

    // int process_event(PHCompositeNode *topNode) override;
  
  protected:
    void cluster(std::vector<RawClusterBuilderHelper::towersStrct> &input_towers, uint caloId) override;
};

#endif
