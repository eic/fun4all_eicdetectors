#ifndef G4CALO__RAWTOWERBUILDERBYHITINDEXBECAL_H
#define G4CALO__RAWTOWERBUILDERBYHITINDEXBECAL_H

#include <calobase/RawTowerDefs.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

/**
 * \brief SubsysReco module creating calorimeter tower objects (RawTowerv1) from hits
 * (PHG4Hit) using j,k indeces of these hits
 *
 */
class RawTowerBuilderByHitIndexBECAL : public SubsysReco
{
 public:
  RawTowerBuilderByHitIndexBECAL(const std::string &name = "RawTowerBuilderByHitIndexBECAL");
  ~RawTowerBuilderByHitIndexBECAL() override {}

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  /** Name of the detector node the G4Hits should be taken from.
   */
  void Detector(const std::string &d);

  /** Specifiy text-file with table for tower mapping
   */
  void GeometryTableFile(const std::string &d)
  {
    m_MappingTowerFile = d;
  }

  /** Define minimum tower energy. After processing an event, towers with lower energy
   * are will be deleted.
   */
  void EminCut(const double e) { m_Emin = e; }

  /** Get prefix for tower collection to identify simulated towers
   * before digitization.
   */
  std::string
  get_sim_tower_node_prefix() const
  {
    return m_SimTowerNodePrefix;
  }

  /** Set prefix for tower collection to identify simulated towers
   * before digitization.
   */
  void
  set_sim_tower_node_prefix(const std::string &simTowerNodePrefix)
  {
    m_SimTowerNodePrefix = simTowerNodePrefix;
  }
  /** Set time window allowed for tower aggregation.
   */
  void
  set_hit_time_window(const double tmax)
  {
    m_Tmax = tmax;
  }
  /** Set gdml geometry loading.
   */
  void
  set_use_gdml(bool usegdml)
  {
    m_useGDML = usegdml;
  }
 private:
  /** Create nodes for output.
   *
   * Name of output node for RawTowerContainer: "TOWER_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);

  /** Read geometry information from table stored in text-file
   */
  bool ReadGeometryFromTable();
  bool ReadGeometryFromGDML();

  RawTowerContainer *m_Towers;
  RawTowerGeomContainer *m_Geoms;

  std::string m_Detector;
  std::string m_SimTowerNodePrefix;

  std::string m_MappingTowerFile;

  RawTowerDefs::CalorimeterId m_CaloId;

  double thickness_wall = -1.;
  double radius = 85.0;
  double tower_length = 45.5;

  double m_Emin;
  double m_Tmax;
  bool m_useGDML;

  std::map<std::string, double> m_GlobalParameterMap;
};

#endif
