#ifndef G4TTL__RAWDIGITBUILDERTTL_H
#define G4TTL__RAWDIGITBUILDERTTL_H

// #include <calobase/RawTowerDefs.h>

#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>               // for G4double, G4int

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
// #include <trackbase/TrkrCluster.h>
#include <map>
#include <string>
#include <utility>

class PHCompositeNode;
// class RawTowerContainer;
// class RawTowerGeomContainer;
// class TrkrHit;
// class TrkrHitSetContainer;
class TrkrClusterContainer;
// class TrkrClusterHitAssoc;
// class TrkrClusterContainer;

/**
 * \brief SubsysReco module creating calorimeter tower objects (RawTowerv1) from hits
 * (PHG4Hit) using j,k indeces of these hits
 *
 */
class RawDigitBuilderTTL : public SubsysReco
{
 public:
  RawDigitBuilderTTL(const std::string &name = "RawDigitBuilderTTL");
  ~RawDigitBuilderTTL() override {}

  //! module initialization
  int Init(PHCompositeNode *topNode) override { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! end of process
  int End(PHCompositeNode *topNode) override { return 0; }

  //! option to turn off z-dimension clustering
  void SetZClustering(const bool make_z_clustering)
  {
    m_makeZClustering = make_z_clustering;
  }
  bool GetZClustering() const
  {
    return m_makeZClustering;
  }

  /** Name of the detector node the G4Hits should be taken from.
   */
  void Detector(const std::string &d);

  /** Define minimum tower energy. After processing an event, towers with lower energy
   * are will be deleted.
   */
  // void EminCut(const double e) { m_Emin = e; }


 private:
  /** Create nodes for output.
   *
   * Name of output node for RawTowerContainer: "TOWER_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);
  void PrintClusters(PHCompositeNode *topNode);
  // void GetPixelGlobalCoordinates(PHG4Hit* g4hit, G4double &xpos, G4double &ypos, G4double &zpos);
  // TrkrHitSetContainer *m_hits;
  TrkrClusterContainer *m_clusterlist; 

  // TrkrClusterHitAssoc *m_clusterhitassoc;

  /** Read geometry information from table stored in text-file
   */

  // RawTowerContainer *m_Towers;
  // RawTowerGeomContainer *m_Geoms;

  std::string m_Detector;
  std::string m_SimTowerNodePrefix;

  // RawTowerDefs::CalorimeterId m_CaloId;

  double m_Emin;

  // settings
  bool m_makeZClustering;  // z_clustering_option
  std::map<std::string, double> m_GlobalParameterMap;
};

#endif
