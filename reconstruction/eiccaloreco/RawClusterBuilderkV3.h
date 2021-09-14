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
    RawClusterBuilderkV3(const std::string &name);

  protected:
    void cluster(std::vector<towersStrct> &input_towers, uint caloId) override;
};

#endif
