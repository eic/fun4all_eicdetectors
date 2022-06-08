#ifndef G4EVAL_EVENTEVALUATOREIC_H
#define G4EVAL_EVENTEVALUATOREIC_H

//===============================================
/// \file EventEvaluatorEIC.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised sPHENIX version)
//===============================================

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>
#include <vector>

class CaloEvalStack;
class PHCompositeNode;
class PHHepMCGenEventMap;
class PHHepMCGenEvent;
class TFile;
class TNtuple;
class TTree;  //Added by Barak

/// \class EventEvaluatorEIC
///
/// \brief Compares reconstructed showers to truth particles
///
/// Plan: This module will trace the reconstructed clusters back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class EventEvaluatorEIC : public SubsysReco
{
 public:
  enum class TrackSource_t : unsigned short
  {
    all = 0,
    inner = 1,
    silicon = 2,
    ttl = 3,
    defecce = 4
  };

  EventEvaluatorEIC(const std::string& name = "EventEvaluatorEIC",
                    const std::string& filename = "g4eval_event.root");
  ~EventEvaluatorEIC() override{};

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_strict(bool b) { _strict = b; }

  void set_do_store_event_level_info(bool b) { _do_store_event_info = b; }
  void set_do_BECAL(bool b) { _do_BECAL = b; }
  void set_do_HCALIN(bool b) { _do_HCALIN = b; }
  void set_do_HCALOUT(bool b) { _do_HCALOUT = b; }
  void set_do_EHCAL(bool b) { _do_EHCAL = b; }
  void set_do_FEMC(bool b) { _do_FEMC = b; }
  void set_do_EEMC(bool b) { _do_EEMC = b; }
  void set_do_DRCALO(bool b) { _do_DRCALO = b; }
  void set_do_LFHCAL(bool b) { _do_LFHCAL = b; }
  void set_do_HITS(bool b) { _do_HITS = b; }
  void set_do_HITS_ABSORBER(bool b) { _do_HITS_ABSORBER = b; }
  void set_do_HITS_CALO(bool b) { _do_HITS_CALO = b; }
  void set_do_TRACKS(bool b) { _do_TRACKS = b; }
  void set_do_PID_LogLikelihood(bool b) { _do_PID_LogLikelihood = b; }
  void set_do_CLUSTERS(bool b) { _do_CLUSTERS = b; }
  void set_do_VERTEX(bool b) { _do_VERTEX = b; }
  void set_do_PROJECTIONS(bool b) { _do_PROJECTIONS = b; }
  void set_do_MCPARTICLES(bool b) { _do_MCPARTICLES = b; }
  void set_do_HEPMC(bool b) { _do_HEPMC = b; }
  void set_do_GEOMETRY(bool b) { _do_GEOMETRY = b; }
  void set_do_BLACKHOLE(bool b) { _do_BLACKHOLE = b; }

  // limit the tracing of towers and clusters back to the truth particles
  // to only those reconstructed objects above a particular energy
  // threshold (evaluation for objects above threshold unaffected)
  void set_reco_tracing_energy_threshold(float thresh, int caloid)
  {
    _reco_e_threshold[caloid] = thresh;
  }
  void set_reco_tracing_energy_thresholdMC(float thresh)
  {
    _reco_e_thresholdMC = thresh;
  }

  //! max depth/generation of the MC_particle/PHG4Particle that would be saved.
  void set_depth_MCstack(int d)
  {
    _depth_MCstack = d;
  }

 private:
  bool _do_store_event_info;
  bool _do_BECAL;
  bool _do_HCALIN;
  bool _do_HCALOUT;
  bool _do_EHCAL;
  bool _do_FEMC;
  bool _do_EEMC;
  bool _do_DRCALO;
  bool _do_LFHCAL;
  bool _do_HITS;
  bool _do_HITS_ABSORBER;
  bool _do_HITS_CALO;
  bool _do_TRACKS;
  bool _do_CLUSTERS;
  bool _do_VERTEX;
  bool _do_PROJECTIONS;
  bool _do_PID_LogLikelihood = false;
  bool _do_MCPARTICLES;
  bool _do_HEPMC;
  bool _do_GEOMETRY;
  bool _do_BLACKHOLE;
  unsigned int _ievent;

  // Event level info
  float _cross_section;
  float _event_weight;
  int _n_generator_accepted;

  // track hits
  int _nHitsLayers;
  int* _hits_layerID;
  int* _hits_trueID;
  float* _hits_x;
  float* _hits_y;
  float* _hits_z;
  float* _hits_x2;
  float* _hits_y2;
  float* _hits_z2;
  float* _hits_t;
  float* _hits_edep;
  float* _hits_lightyield;
  int* _hits_isAbsorber;

  // towers
  int _nTowers_BECAL;
  float* _tower_BECAL_E;
  int* _tower_BECAL_iEta;
  int* _tower_BECAL_iPhi;
  int* _tower_BECAL_trueID;

  // towers
  int _nTowers_HCALIN;
  float* _tower_HCALIN_E;
  int* _tower_HCALIN_iEta;
  int* _tower_HCALIN_iPhi;
  int* _tower_HCALIN_trueID;

  // towers
  int _nTowers_HCALOUT;
  float* _tower_HCALOUT_E;
  int* _tower_HCALOUT_iEta;
  int* _tower_HCALOUT_iPhi;
  int* _tower_HCALOUT_trueID;

  int _nTowers_EHCAL;
  float* _tower_EHCAL_E;
  int* _tower_EHCAL_iEta;
  int* _tower_EHCAL_iPhi;
  int* _tower_EHCAL_trueID;

  int _nTowers_DRCALO;
  float* _tower_DRCALO_E;
  int* _tower_DRCALO_NScint;
  int* _tower_DRCALO_NCerenkov;
  int* _tower_DRCALO_iEta;
  int* _tower_DRCALO_iPhi;
  int* _tower_DRCALO_trueID;

  int _nTowers_LFHCAL;
  float* _tower_LFHCAL_E;
  int* _tower_LFHCAL_iEta;
  int* _tower_LFHCAL_iPhi;
  int* _tower_LFHCAL_iL;
  int* _tower_LFHCAL_trueID;

  int _nTowers_FEMC;
  float* _tower_FEMC_E;
  int* _tower_FEMC_iEta;
  int* _tower_FEMC_iPhi;
  int* _tower_FEMC_trueID;

  int _nTowers_EEMC;
  float* _tower_EEMC_E;
  int* _tower_EEMC_iEta;
  int* _tower_EEMC_iPhi;
  int* _tower_EEMC_trueID;
 
  int _nclusters_HCALIN;
  float* _cluster_HCALIN_E;
  float* _cluster_HCALIN_Eta;
  float* _cluster_HCALIN_Phi;
  int* _cluster_HCALIN_NTower;
  int* _cluster_HCALIN_trueID;

  int _nclusters_HCALOUT;
  float* _cluster_HCALOUT_E;
  float* _cluster_HCALOUT_Eta;
  float* _cluster_HCALOUT_Phi;
  int* _cluster_HCALOUT_NTower;
  int* _cluster_HCALOUT_trueID;

  int _nclusters_EHCAL;
  float* _cluster_EHCAL_E;
  float* _cluster_EHCAL_Eta;
  float* _cluster_EHCAL_Phi;
  int* _cluster_EHCAL_NTower;
  int* _cluster_EHCAL_trueID;

  int _nclusters_FEMC;
  float* _cluster_FEMC_E;
  float* _cluster_FEMC_Eta;
  float* _cluster_FEMC_Phi;
  int* _cluster_FEMC_NTower;
  int* _cluster_FEMC_trueID;

  int _nclusters_EEMC;
  float* _cluster_EEMC_E;
  float* _cluster_EEMC_Eta;
  float* _cluster_EEMC_Phi;
  int* _cluster_EEMC_NTower;
  int* _cluster_EEMC_trueID;

  // vertex
  float _vertex_x;
  float _vertex_y;
  float _vertex_z;
  int _vertex_NCont;
  float _vertex_true_x;
  float _vertex_true_y;
  float _vertex_true_z;

  // tracks
  int _nTracks;
  float* _track_ID;
  short* _track_charge;
  float* _track_px;
  float* _track_py;
  float* _track_pz;
  float* _track_x;
  float* _track_y;
  float* _track_z;
  float* _track_ndf;
  short* _track_nHits;
  unsigned int* _track_hitsEncoded;
  float* _track_chi2;
  float* _track_dca;
  float* _track_dca_2d;
  float* _track_trueID;
  unsigned short* _track_source;

  // log likelihood summary for PID detectors, per track information
  std::vector<float> _track_pion_LL;
  std::vector<float> _track_kaon_LL;
  std::vector<float> _track_proton_LL;

  int _nProjections;
  float* _track_ProjTrackID;
  int* _track_ProjLayer;
  float* _track_TLP_x;
  float* _track_TLP_y;
  float* _track_TLP_z;
  float* _track_TLP_t;
  float* _track_TLP_px;
  float* _track_TLP_py;
  float* _track_TLP_pz;
  float* _track_TLP_true_x;
  float* _track_TLP_true_y;
  float* _track_TLP_true_z;
  float* _track_TLP_true_t;

  // MC particles
  int _nMCPart;
  int* _mcpart_ID;
  int* _mcpart_ID_parent;
  int* _mcpart_PDG;
  float* _mcpart_E;
  float* _mcpart_px;
  float* _mcpart_py;
  float* _mcpart_pz;
  float* _mcpart_x;
  float* _mcpart_y;
  float* _mcpart_z;
  int* _mcpart_BCID;

  // MC particles
  int _nHepmcp;
  int _hepmcp_procid;
  float _hepmcp_x1;
  float _hepmcp_x2;
  float _hepmcp_Q2;
  float _hepmcp_vtx_x;
  float _hepmcp_vtx_y;
  float _hepmcp_vtx_z;
  float _hepmcp_vtx_t;
  //  float* _hepmcp_ID_parent;
  int* _hepmcp_status;
  int* _hepmcp_PDG;
  float* _hepmcp_E;
  float* _hepmcp_px;
  float* _hepmcp_py;
  float* _hepmcp_pz;
  int* _hepmcp_m1;
  int* _hepmcp_m2;
  int* _hepmcp_BCID;

  int _calo_ID;
  int _calo_towers_N;
  int* _calo_towers_iEta;
  int* _calo_towers_iPhi;
  int* _calo_towers_iL;
  float* _calo_towers_Eta;
  float* _calo_towers_Phi;
  float* _calo_towers_x;
  float* _calo_towers_y;
  float* _calo_towers_z;
  int* _geometry_done;

  float* _reco_e_threshold;
  float _reco_e_thresholdMC;
  int _depth_MCstack;

  CaloEvalStack* _caloevalstackBECAL;
  CaloEvalStack* _caloevalstackHCALIN;
  CaloEvalStack* _caloevalstackHCALOUT;
  CaloEvalStack* _caloevalstackEHCAL;
  CaloEvalStack* _caloevalstackDRCALO;
  CaloEvalStack* _caloevalstackLFHCAL;
  CaloEvalStack* _caloevalstackFEMC;
  CaloEvalStack* _caloevalstackEEMC;

  //----------------------------------
  // evaluator output ntuples

  bool _strict;

  TTree* _event_tree;     //Added by Barak
  TTree* _geometry_tree;  //Added by Barak

  // evaluator output file
  std::string _filename;
  TFile* _tfile;
  TFile* _tfile_geometry;

  // subroutines
  int GetProjectionIndex(std::string projname);           ///< return track projection index for given track projection layer
  std::string GetProjectionNameFromIndex(int projindex);  ///< return track projection layer name from projection index (see GetProjectionIndex)
  int GetExponentFromProjectionIndex(int projindex);      ///< return exponent for bitwise encoding of tracking layers
  void fillOutputNtuples(PHCompositeNode* topNode);       ///< dump the evaluator information into ntuple for external analysis
  void resetGeometryArrays();                             ///< reset the tree variables before filling for a new event
  void resetBuffer();                                     ///< reset the tree variables before filling for a new event

  const int _maxNHits = 1000000;
  const int _maxNTowers = 50 * 50;
  const int _maxNTowersCentral = 2000;
  const int _maxNTowersDR = 3000 * 3000;
  const int _maxNTowersCalo = 5000000;
  const int _maxNclusters = 100;
  const int _maxNclustersCentral = 2000;
  const int _maxNTracks = 200;
  const int _maxNProjections = 2000;
  const int _maxNMCPart = 100000;
  const int _maxNHepmcp = 1000;
  const int _maxNCalo = 15;

  enum calotype
  {
    kFEMC = 1,
    kDRCALO = 2,
    kEEMC = 3,
    kEHCAL = 5,
    kHCALIN = 6,
    kHCALOUT = 7,
    kLFHCAL = 8,
    kBECAL = 10,
  };
};

#endif  // G4EVAL_EVENTEVALUATOR_H
