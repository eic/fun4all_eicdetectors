#ifndef EICCALORECO_RAWCLUSTERBUILDERKMA_H
#define EICCALORECO_RAWCLUSTERBUILDERKMA_H

#include <RawClusterBuilderHelper.h>

#include <string>
#include <vector>

class RawClusterBuilderkMA : public RawClusterBuilderHelper {
  public:
    RawClusterBuilderkMA(const std::string &name);

  protected:
    void cluster(std::vector<RawClusterBuilderHelper::towersStrct> &input_towers, uint caloId) override;
};

#endif