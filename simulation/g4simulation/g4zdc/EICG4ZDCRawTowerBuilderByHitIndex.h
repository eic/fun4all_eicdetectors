#ifndef EICG4ZDC_RAWTOWERBUILDER_H
#define EICG4ZDC_RAWTOWERBUILDER_H

#include <eiczdcbase/RawTowerZDCDefs.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class PHCompositeNode;
class RawTowerZDCContainer;
class RawTowerZDCGeomContainer;

/**
 * \brief SubsysReco module creating calorimeter tower objects (RawTowerZDCv1) from hits
 * (PHG4Hit) using j,k,l indeces of these hits
 *
 */
class EICG4ZDCRawTowerBuilderByHitIndex : public SubsysReco
{
 public:
  EICG4ZDCRawTowerBuilderByHitIndex(const std::string &name = "EICG4ZDCRawTowerBuilderByHitIndex");
  ~EICG4ZDCRawTowerBuilderByHitIndex() override {}

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  /** Name of the detector node the G4Hits should be taken from.
   */
  void Detector(const std::string &d);

  void SubDetector(const std::string &d);
  
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

 private:
  /** Create nodes for output.
   *
   * Name of output node for RawTowerZDCContainer: "TOWER_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);

  /** Read geometry information from table stored in text-file
   */
  bool ReadGeometryFromTable();

  RawTowerZDCContainer *m_Towers;
  RawTowerZDCGeomContainer *m_Geoms;

  std::string m_Detector;
  std::string m_SubDetector;
  std::string m_SimTowerNodePrefix;

  std::string m_MappingTowerFile;

  int m_SubDetID;
 
  RawTowerZDCDefs::CalorimeterId m_CaloId;

  double m_Emin;
  double m_TowerDepth;
  double m_ThicknessAbsorber;
  double  m_ThicknessScintilator;
  std::map<std::string, double> m_GlobalParameterMap;
  std::map<int, int> m_TowerIDtoLayerIDMap;
};

#endif
