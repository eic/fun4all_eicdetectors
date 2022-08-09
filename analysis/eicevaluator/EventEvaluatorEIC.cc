#include "EventEvaluatorEIC.h"

#include "g4eval/CaloEvalStack.h"
#include "g4eval/CaloRawClusterEval.h"
#include "g4eval/CaloRawTowerEval.h"
#include "g4eval/CaloTruthEval.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxVertex.h>  // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv2.h>

#include <eicpidbase/EICPIDParticle.h>
#include <eicpidbase/EICPIDParticleContainer.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHGenIntegral.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <utility>
// #include <fstream>

using namespace std;

EventEvaluatorEIC::EventEvaluatorEIC(const string& name, const string& filename)
  : SubsysReco(name)
  , _do_store_event_info(false)
  , _do_BECAL(false)
  , _do_HCALIN(false)
  , _do_HCALOUT(false)
  , _do_EHCAL(false)
  , _do_FEMC(false)
  , _do_EEMC(false)
  , _do_DRCALO(false)
  , _do_LFHCAL(false)
  , _do_HITS(false)
  , _do_HITS_ABSORBER(false)
  , _do_HITS_CALO(false)
  , _do_TRACKS(false)
  , _do_CLUSTERS(false)
  , _do_VERTEX(false)
  , _do_PROJECTIONS(false)
  , _do_MCPARTICLES(false)
  , _do_HEPMC(false)
  , _do_GEOMETRY(false)
  , _do_BLACKHOLE(false)
  , _ievent(0)
  , _cross_section(0)
  , _event_weight(0)
  , _n_generator_accepted(0)
  , _nHitsLayers(0)
  , _hits_layerID(0)
  , _hits_trueID(0)
  , _hits_x(0)
  , _hits_y(0)
  , _hits_z(0)
  , _hits_x2(0)
  , _hits_y2(0)
  , _hits_z2(0)
  , _hits_t(0)
  , _hits_edep(0)
  , _hits_lightyield(0)
  , _hits_isAbsorber(0)

  , _nTowers_BECAL(0)
  , _tower_BECAL_E(0)
  , _tower_BECAL_iEta(0)
  , _tower_BECAL_iPhi(0)
  , _tower_BECAL_trueID(0)

  , _nTowers_HCALIN(0)
  , _tower_HCALIN_E(0)
  , _tower_HCALIN_iEta(0)
  , _tower_HCALIN_iPhi(0)
  , _tower_HCALIN_trueID(0)

  , _nTowers_HCALOUT(0)
  , _tower_HCALOUT_E(0)
  , _tower_HCALOUT_iEta(0)
  , _tower_HCALOUT_iPhi(0)
  , _tower_HCALOUT_trueID(0)

  , _nTowers_EHCAL(0)
  , _tower_EHCAL_E(0)
  , _tower_EHCAL_iEta(0)
  , _tower_EHCAL_iPhi(0)
  , _tower_EHCAL_trueID(0)

  , _nTowers_DRCALO(0)
  , _tower_DRCALO_E(0)
  , _tower_DRCALO_NScint(0)
  , _tower_DRCALO_NCerenkov(0)
  , _tower_DRCALO_iEta(0)
  , _tower_DRCALO_iPhi(0)
  , _tower_DRCALO_trueID(0)

  , _nTowers_LFHCAL(0)
  , _tower_LFHCAL_E(0)
  , _tower_LFHCAL_iEta(0)
  , _tower_LFHCAL_iPhi(0)
  , _tower_LFHCAL_iL(0)
  , _tower_LFHCAL_trueID(0)

  , _nTowers_FEMC(0)
  , _tower_FEMC_E(0)
  , _tower_FEMC_iEta(0)
  , _tower_FEMC_iPhi(0)
  , _tower_FEMC_trueID(0)

  , _nTowers_EEMC(0)
  , _tower_EEMC_E(0)
  , _tower_EEMC_iEta(0)
  , _tower_EEMC_iPhi(0)
  , _tower_EEMC_trueID(0)

  , _nclusters_HCALIN(0)
  , _cluster_HCALIN_E(0)
  , _cluster_HCALIN_Eta(0)
  , _cluster_HCALIN_Phi(0)
  , _cluster_HCALIN_NTower(0)
  , _cluster_HCALIN_trueID(0)

  , _nclusters_HCALOUT(0)
  , _cluster_HCALOUT_E(0)
  , _cluster_HCALOUT_Eta(0)
  , _cluster_HCALOUT_Phi(0)
  , _cluster_HCALOUT_NTower(0)
  , _cluster_HCALOUT_trueID(0)

  , _nclusters_EHCAL(0)
  , _cluster_EHCAL_E(0)
  , _cluster_EHCAL_Eta(0)
  , _cluster_EHCAL_Phi(0)
  , _cluster_EHCAL_NTower(0)
  , _cluster_EHCAL_trueID(0)

  , _nclusters_FEMC(0)
  , _cluster_FEMC_E(0)
  , _cluster_FEMC_Eta(0)
  , _cluster_FEMC_Phi(0)
  , _cluster_FEMC_NTower(0)
  , _cluster_FEMC_trueID(0)

  , _nclusters_EEMC(0)
  , _cluster_EEMC_E(0)
  , _cluster_EEMC_Eta(0)
  , _cluster_EEMC_Phi(0)
  , _cluster_EEMC_NTower(0)
  , _cluster_EEMC_trueID(0)

  , _vertex_x(0)
  , _vertex_y(0)
  , _vertex_z(0)
  , _vertex_NCont(0)
  , _vertex_true_x(0)
  , _vertex_true_y(0)
  , _vertex_true_z(0)

  , _nTracks(0)
  , _track_ID(0)
  , _track_charge(0)
  , _track_px(0)
  , _track_py(0)
  , _track_pz(0)
  , _track_x(0)
  , _track_y(0)
  , _track_z(0)
  , _track_ndf(0)
  , _track_nHits(0)
  , _track_hitsEncoded(0)
  , _track_chi2(0)
  , _track_dca(0)
  , _track_dca_2d(0)
  , _track_trueID(0)
  , _track_source(0)
  , _nProjections(0)
  , _track_ProjTrackID(0)
  , _track_ProjLayer(0)
  , _track_TLP_x(0)
  , _track_TLP_y(0)
  , _track_TLP_z(0)
  , _track_TLP_t(0)
  , _track_TLP_px(0)
  , _track_TLP_py(0)
  , _track_TLP_pz(0)
  , _track_TLP_true_x(0)
  , _track_TLP_true_y(0)
  , _track_TLP_true_z(0)
  , _track_TLP_true_t(0)

  , _nMCPart(0)
  , _mcpart_ID(0)
  , _mcpart_ID_parent(0)
  , _mcpart_PDG(0)
  , _mcpart_E(0)
  , _mcpart_px(0)
  , _mcpart_py(0)
  , _mcpart_pz(0)
  , _mcpart_x(0)
  , _mcpart_y(0)
  , _mcpart_z(0)
  , _mcpart_BCID(0)

  , _nHepmcp(0)
  , _hepmcp_procid(0)
  , _hepmcp_x1(NAN)
  , _hepmcp_x2(NAN)
  , _hepmcp_Q2(NAN)
  , _hepmcp_vtx_x(NAN)
  , _hepmcp_vtx_y(NAN)
  , _hepmcp_vtx_z(NAN)
  , _hepmcp_vtx_t(NAN)

  //  , _hepmcp_ID_parent(0)
  , _hepmcp_status(0)
  , _hepmcp_PDG(0)
  , _hepmcp_E(0)
  , _hepmcp_px(0)
  , _hepmcp_py(0)
  , _hepmcp_pz(0)
  , _hepmcp_m1(0)
  , _hepmcp_m2(0)
  , _hepmcp_BCID(0)

  , _calo_ID(0)
  , _calo_towers_N(0)
  , _calo_towers_iEta(0)
  , _calo_towers_iPhi(0)
  , _calo_towers_iL(0)
  , _calo_towers_Eta(0)
  , _calo_towers_Phi(0)
  , _calo_towers_x(0)
  , _calo_towers_y(0)
  , _calo_towers_z(0)
  , _geometry_done(0)

  , _reco_e_threshold(0)
  , _reco_e_thresholdMC(0.001)
  , _depth_MCstack(0)
  , _caloevalstackBECAL(nullptr)
  , _caloevalstackHCALIN(nullptr)
  , _caloevalstackHCALOUT(nullptr)
  , _caloevalstackEHCAL(nullptr)
  , _caloevalstackDRCALO(nullptr)
  , _caloevalstackLFHCAL(nullptr)
  , _caloevalstackFEMC(nullptr)
  , _caloevalstackEEMC(nullptr)
  , _strict(false)
  , _event_tree(nullptr)
  , _geometry_tree(nullptr)
  , _filename(filename)
  , _tfile(nullptr)
  , _tfile_geometry(nullptr)
{
  _reco_e_threshold = new float[_maxNCalo];
  _reco_e_threshold[kFEMC] = 0.005;
  _reco_e_threshold[kDRCALO] = 0.0;
  _reco_e_threshold[kEEMC] = 0.005;
  _reco_e_threshold[kEHCAL] = 0.05;
  _reco_e_threshold[kHCALIN] = 0.01;
  _reco_e_threshold[kHCALOUT] = 0.05;
  _reco_e_threshold[kLFHCAL] = 0.001;
  _reco_e_threshold[kBECAL] = 0.001;

  _hits_layerID = new int[_maxNHits];
  _hits_trueID = new int[_maxNHits];
  _hits_x = new float[_maxNHits];
  _hits_y = new float[_maxNHits];
  _hits_z = new float[_maxNHits];
  _hits_x2 = new float[_maxNHits];
  _hits_y2 = new float[_maxNHits];
  _hits_z2 = new float[_maxNHits];
  _hits_t = new float[_maxNHits];
  _hits_edep = new float[_maxNHits];
  _hits_lightyield = new float[_maxNHits];
  _hits_isAbsorber = new int[_maxNHits];

  _tower_BECAL_E = new float[_maxNTowers];
  _tower_BECAL_iEta = new int[_maxNTowers];
  _tower_BECAL_iPhi = new int[_maxNTowers];
  _tower_BECAL_trueID = new int[_maxNTowers];

  _tower_HCALIN_E = new float[_maxNTowersCentral];
  _tower_HCALIN_iEta = new int[_maxNTowersCentral];
  _tower_HCALIN_iPhi = new int[_maxNTowersCentral];
  _tower_HCALIN_trueID = new int[_maxNTowersCentral];
  _cluster_HCALIN_E = new float[_maxNclusters];
  _cluster_HCALIN_Eta = new float[_maxNclusters];
  _cluster_HCALIN_Phi = new float[_maxNclusters];
  _cluster_HCALIN_NTower = new int[_maxNclusters];
  _cluster_HCALIN_trueID = new int[_maxNclusters];

  _tower_HCALOUT_E = new float[_maxNTowersCentral];
  _tower_HCALOUT_iEta = new int[_maxNTowersCentral];
  _tower_HCALOUT_iPhi = new int[_maxNTowersCentral];
  _tower_HCALOUT_trueID = new int[_maxNTowersCentral];
  _cluster_HCALOUT_E = new float[_maxNclusters];
  _cluster_HCALOUT_Eta = new float[_maxNclusters];
  _cluster_HCALOUT_Phi = new float[_maxNclusters];
  _cluster_HCALOUT_NTower = new int[_maxNclusters];
  _cluster_HCALOUT_trueID = new int[_maxNclusters];

  _tower_EHCAL_E = new float[_maxNTowers];
  _tower_EHCAL_iEta = new int[_maxNTowers];
  _tower_EHCAL_iPhi = new int[_maxNTowers];
  _tower_EHCAL_trueID = new int[_maxNTowers];
  _cluster_EHCAL_E = new float[_maxNclusters];
  _cluster_EHCAL_Eta = new float[_maxNclusters];
  _cluster_EHCAL_Phi = new float[_maxNclusters];
  _cluster_EHCAL_NTower = new int[_maxNclusters];
  _cluster_EHCAL_trueID = new int[_maxNclusters];

  _tower_DRCALO_E = new float[_maxNTowersDR];
  _tower_DRCALO_NScint = new int[_maxNTowersDR];
  _tower_DRCALO_NCerenkov = new int[_maxNTowersDR];
  _tower_DRCALO_iEta = new int[_maxNTowersDR];
  _tower_DRCALO_iPhi = new int[_maxNTowersDR];
  _tower_DRCALO_trueID = new int[_maxNTowersDR];

  _tower_LFHCAL_E = new float[_maxNTowers];
  _tower_LFHCAL_iEta = new int[_maxNTowers];
  _tower_LFHCAL_iPhi = new int[_maxNTowers];
  _tower_LFHCAL_iL = new int[_maxNTowers];
  _tower_LFHCAL_trueID = new int[_maxNTowers];

  _tower_FEMC_E = new float[_maxNTowers];
  _tower_FEMC_iEta = new int[_maxNTowers];
  _tower_FEMC_iPhi = new int[_maxNTowers];
  _tower_FEMC_trueID = new int[_maxNTowers];
  _cluster_FEMC_E = new float[_maxNclusters];
  _cluster_FEMC_Eta = new float[_maxNclusters];
  _cluster_FEMC_Phi = new float[_maxNclusters];
  _cluster_FEMC_NTower = new int[_maxNclusters];
  _cluster_FEMC_trueID = new int[_maxNclusters];

  _tower_EEMC_E = new float[_maxNTowers];
  _tower_EEMC_iEta = new int[_maxNTowers];
  _tower_EEMC_iPhi = new int[_maxNTowers];
  _tower_EEMC_trueID = new int[_maxNTowers];
  _cluster_EEMC_E = new float[_maxNclusters];
  _cluster_EEMC_Eta = new float[_maxNclusters];
  _cluster_EEMC_Phi = new float[_maxNclusters];
  _cluster_EEMC_NTower = new int[_maxNclusters];
  _cluster_EEMC_trueID = new int[_maxNclusters];

  _track_ID = new float[_maxNTracks];
  _track_charge = new short[_maxNTracks];
  _track_trueID = new float[_maxNTracks];
  _track_px = new float[_maxNTracks];
  _track_py = new float[_maxNTracks];
  _track_pz = new float[_maxNTracks];
  _track_x = new float[_maxNTracks];
  _track_y = new float[_maxNTracks];
  _track_z = new float[_maxNTracks];
  _track_ndf = new float[_maxNTracks];
  _track_nHits = new short[_maxNTracks];
  _track_hitsEncoded = new unsigned int[_maxNTracks];
  _track_chi2 = new float[_maxNTracks];
  _track_dca = new float[_maxNTracks];
  _track_dca_2d = new float[_maxNTracks];
  _track_source = new unsigned short[_maxNTracks];
  _track_ProjTrackID = new float[_maxNProjections];
  _track_ProjLayer = new int[_maxNProjections];
  _track_TLP_x = new float[_maxNProjections];
  _track_TLP_y = new float[_maxNProjections];
  _track_TLP_z = new float[_maxNProjections];
  _track_TLP_t = new float[_maxNProjections];
  _track_TLP_px = new float[_maxNProjections];
  _track_TLP_py = new float[_maxNProjections];
  _track_TLP_pz = new float[_maxNProjections];
  _track_TLP_true_x = new float[_maxNProjections];
  _track_TLP_true_y = new float[_maxNProjections];
  _track_TLP_true_z = new float[_maxNProjections];
  _track_TLP_true_t = new float[_maxNProjections];
  _track_pion_LL.resize(_maxNTracks, -100);
  _track_kaon_LL.resize(_maxNTracks, -100);
  _track_proton_LL.resize(_maxNTracks, -100);

  _mcpart_ID = new int[_maxNMCPart];
  _mcpart_ID_parent = new int[_maxNMCPart];
  _mcpart_PDG = new int[_maxNMCPart];
  _mcpart_E = new float[_maxNMCPart];
  _mcpart_px = new float[_maxNMCPart];
  _mcpart_py = new float[_maxNMCPart];
  _mcpart_pz = new float[_maxNMCPart];
  _mcpart_x = new float[_maxNMCPart];
  _mcpart_y = new float[_maxNMCPart];
  _mcpart_z = new float[_maxNMCPart];
  _mcpart_BCID = new int[_maxNMCPart];

  _hepmcp_BCID = new int[_maxNHepmcp];
  //  _hepmcp_ID_parent = new float[_maxNHepmcp];
  _hepmcp_status = new int[_maxNHepmcp];
  _hepmcp_PDG = new int[_maxNHepmcp];
  _hepmcp_E = new float[_maxNHepmcp];
  _hepmcp_px = new float[_maxNHepmcp];
  _hepmcp_py = new float[_maxNHepmcp];
  _hepmcp_pz = new float[_maxNHepmcp];
  _hepmcp_m1 = new int[_maxNHepmcp];
  _hepmcp_m2 = new int[_maxNHepmcp];

  _calo_towers_iEta = new int[_maxNTowersCalo];
  _calo_towers_iPhi = new int[_maxNTowersCalo];
  _calo_towers_iL = new int[_maxNTowersCalo];
  _calo_towers_Eta = new float[_maxNTowersCalo];
  _calo_towers_Phi = new float[_maxNTowersCalo];
  _calo_towers_x = new float[_maxNTowersCalo];
  _calo_towers_y = new float[_maxNTowersCalo];
  _calo_towers_z = new float[_maxNTowersCalo];
  _geometry_done = new int[20];
  for (int igem = 0; igem < 20; igem++) _geometry_done[igem] = 0;
}

int EventEvaluatorEIC::Init(PHCompositeNode* topNode)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

  _event_tree = new TTree("event_tree", "event_tree");
  if (_do_store_event_info)
  {
    // Event level info. This isn't the most efficient way to store this info, but it's straightforward
    // within the structure of the class, so the size is small compared to the rest of the output.
    _event_tree->Branch("cross_section", &_cross_section, "cross_section/F");
    _event_tree->Branch("event_weight", &_event_weight, "event_weight/F");
    _event_tree->Branch("n_generator_accepted", &_n_generator_accepted, "n_generator_accepted/I");
  }
  // tracks and hits
  if (_do_HITS)
  {
    _event_tree->Branch("nHits", &_nHitsLayers, "nHits/I");
    _event_tree->Branch("hits_layerID", _hits_layerID, "hits_layerID[nHits]/I");
    _event_tree->Branch("hits_trueID", _hits_trueID, "hits_trueID[nHits]/I");
    _event_tree->Branch("hits_x", _hits_x, "hits_x[nHits]/F");
    _event_tree->Branch("hits_y", _hits_y, "hits_y[nHits]/F");
    _event_tree->Branch("hits_z", _hits_z, "hits_z[nHits]/F");
    _event_tree->Branch("hits_x2", _hits_x2, "hits_x2[nHits]/F");
    _event_tree->Branch("hits_y2", _hits_y2, "hits_y2[nHits]/F");
    _event_tree->Branch("hits_z2", _hits_z2, "hits_z2[nHits]/F");
    _event_tree->Branch("hits_t", _hits_t, "hits_t[nHits]/F");
    _event_tree->Branch("hits_edep", _hits_edep, "hits_edep[nHits]/F");
    _event_tree->Branch("hits_lightyield", _hits_lightyield, "hits_lightyield[nHits]/F");
    _event_tree->Branch("hits_isAbsorber", _hits_isAbsorber, "hits_isAbsorber[nHits]/I");
  }

  if (_do_TRACKS)
  {
    _event_tree->Branch("nTracks", &_nTracks, "nTracks/I");
    _event_tree->Branch("tracks_ID", _track_ID, "tracks_ID[nTracks]/F");
    _event_tree->Branch("tracks_charge", _track_charge, "tracks_charge[nTracks]/S");
    _event_tree->Branch("tracks_px", _track_px, "tracks_px[nTracks]/F");
    _event_tree->Branch("tracks_py", _track_py, "tracks_py[nTracks]/F");
    _event_tree->Branch("tracks_pz", _track_pz, "tracks_pz[nTracks]/F");
    _event_tree->Branch("tracks_x", _track_x, "tracks_x[nTracks]/F");
    _event_tree->Branch("tracks_y", _track_y, "tracks_y[nTracks]/F");
    _event_tree->Branch("tracks_z", _track_z, "tracks_z[nTracks]/F");
    _event_tree->Branch("tracks_ndf", _track_ndf, "tracks_ndf[nTracks]/F");
    _event_tree->Branch("tracks_chi2", _track_chi2, "tracks_chi2[nTracks]/F");
    _event_tree->Branch("tracks_nHits", _track_nHits, "tracks_nHits[nTracks]/S");
    _event_tree->Branch("tracks_hitsEncoded", _track_hitsEncoded, "tracks_hitsEncoded[nTracks]/i");
    _event_tree->Branch("tracks_dca", _track_dca, "tracks_dca[nTracks]/F");
    _event_tree->Branch("tracks_dca_2d", _track_dca_2d, "tracks_dca_2d[nTracks]/F");
    _event_tree->Branch("tracks_trueID", _track_trueID, "tracks_trueID[nTracks]/F");
    _event_tree->Branch("tracks_source", _track_source, "tracks_source[nTracks]/s");

    if (_do_PID_LogLikelihood)
    {
      // save hadron PID log likelihood
      _event_tree->Branch("track_pion_LL", _track_pion_LL.data(), "track_pion_LL[nTracks]/F");
      _event_tree->Branch("track_kaon_LL", _track_kaon_LL.data(), "track_kaon_LL[nTracks]/F");
      _event_tree->Branch("track_proton_LL", _track_proton_LL.data(), "track_proton_LL[nTracks]/F");
    }
  }
  if (_do_PROJECTIONS)
  {
    _event_tree->Branch("nProjections", &_nProjections, "nProjections/I");
    _event_tree->Branch("track_ProjTrackID", _track_ProjTrackID, "track_ProjTrackID[nProjections]/F");
    _event_tree->Branch("track_ProjLayer", _track_ProjLayer, "track_ProjLayer[nProjections]/I");
    _event_tree->Branch("track_TLP_x", _track_TLP_x, "track_TLP_x[nProjections]/F");
    _event_tree->Branch("track_TLP_y", _track_TLP_y, "track_TLP_y[nProjections]/F");
    _event_tree->Branch("track_TLP_z", _track_TLP_z, "track_TLP_z[nProjections]/F");
    _event_tree->Branch("track_TLP_t", _track_TLP_t, "track_TLP_t[nProjections]/F");
    _event_tree->Branch("track_TLP_px", _track_TLP_px, "track_TLP_px[nProjections]/F");
    _event_tree->Branch("track_TLP_py", _track_TLP_py, "track_TLP_py[nProjections]/F");
    _event_tree->Branch("track_TLP_pz", _track_TLP_pz, "track_TLP_pz[nProjections]/F");
    _event_tree->Branch("track_TLP_true_x", _track_TLP_true_x, "track_TLP_true_x[nProjections]/F");
    _event_tree->Branch("track_TLP_true_y", _track_TLP_true_y, "track_TLP_true_y[nProjections]/F");
    _event_tree->Branch("track_TLP_true_z", _track_TLP_true_z, "track_TLP_true_z[nProjections]/F");
    _event_tree->Branch("track_TLP_true_t", _track_TLP_true_t, "track_TLP_true_t[nProjections]/F");
  }
  if (_do_BECAL)
  {
    // towers BECAL
    _event_tree->Branch("tower_BECAL_N", &_nTowers_BECAL, "tower_BECAL_N/I");
    _event_tree->Branch("tower_BECAL_E", _tower_BECAL_E, "tower_BECAL_E[tower_BECAL_N]/F");
    _event_tree->Branch("tower_BECAL_iEta", _tower_BECAL_iEta, "tower_BECAL_iEta[tower_BECAL_N]/I");
    _event_tree->Branch("tower_BECAL_iPhi", _tower_BECAL_iPhi, "tower_BECAL_iPhi[tower_BECAL_N]/I");
    _event_tree->Branch("tower_BECAL_trueID", _tower_BECAL_trueID, "tower_BECAL_trueID[tower_BECAL_N]/I");
  }
  if (_do_HCALIN)
  {
    // towers HCAL-in
    _event_tree->Branch("tower_HCALIN_N", &_nTowers_HCALIN, "tower_HCALIN_N/I");
    _event_tree->Branch("tower_HCALIN_E", _tower_HCALIN_E, "tower_HCALIN_E[tower_HCALIN_N]/F");
    _event_tree->Branch("tower_HCALIN_iEta", _tower_HCALIN_iEta, "tower_HCALIN_iEta[tower_HCALIN_N]/I");
    _event_tree->Branch("tower_HCALIN_iPhi", _tower_HCALIN_iPhi, "tower_HCALIN_iPhi[tower_HCALIN_N]/I");
    _event_tree->Branch("tower_HCALIN_trueID", _tower_HCALIN_trueID, "tower_HCALIN_trueID[tower_HCALIN_N]/I");
    if (_do_CLUSTERS)
    {
      // clusters HCAL-in
      _event_tree->Branch("cluster_HCALIN_N", &_nclusters_HCALIN, "cluster_HCALIN_N/I");
      _event_tree->Branch("cluster_HCALIN_E", _cluster_HCALIN_E, "cluster_HCALIN_E[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_Eta", _cluster_HCALIN_Eta, "cluster_HCALIN_Eta[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_Phi", _cluster_HCALIN_Phi, "cluster_HCALIN_Phi[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_NTower", _cluster_HCALIN_NTower, "cluster_HCALIN_NTower[cluster_HCALIN_N]/I");
      _event_tree->Branch("cluster_HCALIN_trueID", _cluster_HCALIN_trueID, "cluster_HCALIN_trueID[cluster_HCALIN_N]/I");
    }
  }
  if (_do_HCALOUT)
  {
    // towers HCAL-out
    _event_tree->Branch("tower_HCALOUT_N", &_nTowers_HCALOUT, "tower_HCALOUT_N/I");
    _event_tree->Branch("tower_HCALOUT_E", _tower_HCALOUT_E, "tower_HCALOUT_E[tower_HCALOUT_N]/F");
    _event_tree->Branch("tower_HCALOUT_iEta", _tower_HCALOUT_iEta, "tower_HCALOUT_iEta[tower_HCALOUT_N]/I");
    _event_tree->Branch("tower_HCALOUT_iPhi", _tower_HCALOUT_iPhi, "tower_HCALOUT_iPhi[tower_HCALOUT_N]/I");
    _event_tree->Branch("tower_HCALOUT_trueID", _tower_HCALOUT_trueID, "tower_HCALOUT_trueID[tower_HCALOUT_N]/I");
    if (_do_CLUSTERS)
    {
      // clusters HCAL-out
      _event_tree->Branch("cluster_HCALOUT_N", &_nclusters_HCALOUT, "cluster_HCALOUT_N/I");
      _event_tree->Branch("cluster_HCALOUT_E", _cluster_HCALOUT_E, "cluster_HCALOUT_E[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_Eta", _cluster_HCALOUT_Eta, "cluster_HCALOUT_Eta[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_Phi", _cluster_HCALOUT_Phi, "cluster_HCALOUT_Phi[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_NTower", _cluster_HCALOUT_NTower, "cluster_HCALOUT_NTower[cluster_HCALOUT_N]/I");
      _event_tree->Branch("cluster_HCALOUT_trueID", _cluster_HCALOUT_trueID, "cluster_HCALOUT_trueID[cluster_HCALOUT_N]/I");
    }
  }
  if (_do_EHCAL)
  {
    // towers EHCAL
    _event_tree->Branch("tower_EHCAL_N", &_nTowers_EHCAL, "tower_EHCAL_N/I");
    _event_tree->Branch("tower_EHCAL_E", _tower_EHCAL_E, "tower_EHCAL_E[tower_EHCAL_N]/F");
    _event_tree->Branch("tower_EHCAL_iEta", _tower_EHCAL_iEta, "tower_EHCAL_iEta[tower_EHCAL_N]/I");
    _event_tree->Branch("tower_EHCAL_iPhi", _tower_EHCAL_iPhi, "tower_EHCAL_iPhi[tower_EHCAL_N]/I");
    _event_tree->Branch("tower_EHCAL_trueID", _tower_EHCAL_trueID, "tower_EHCAL_trueID[tower_EHCAL_N]/I");
    if (_do_CLUSTERS)
    {
      // clusters EHCAL
      _event_tree->Branch("cluster_EHCAL_N", &_nclusters_EHCAL, "cluster_EHCAL_N/I");
      _event_tree->Branch("cluster_EHCAL_E", _cluster_EHCAL_E, "cluster_EHCAL_E[cluster_EHCAL_N]/F");
      _event_tree->Branch("cluster_EHCAL_Eta", _cluster_EHCAL_Eta, "cluster_EHCAL_Eta[cluster_EHCAL_N]/F");
      _event_tree->Branch("cluster_EHCAL_Phi", _cluster_EHCAL_Phi, "cluster_EHCAL_Phi[cluster_EHCAL_N]/F");
      _event_tree->Branch("cluster_EHCAL_NTower", _cluster_EHCAL_NTower, "cluster_EHCAL_NTower[cluster_EHCAL_N]/I");
      _event_tree->Branch("cluster_EHCAL_trueID", _cluster_EHCAL_trueID, "cluster_EHCAL_trueID[cluster_EHCAL_N]/I");
    }
  }
  if (_do_DRCALO)
  {
    // towers DRCALO
    _event_tree->Branch("tower_DRCALO_N", &_nTowers_DRCALO, "tower_DRCALO_N/I");
    _event_tree->Branch("tower_DRCALO_E", _tower_DRCALO_E, "tower_DRCALO_E[tower_DRCALO_N]/F");
    _event_tree->Branch("tower_DRCALO_NScint", _tower_DRCALO_NScint, "tower_DRCALO_NScint[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_NCerenkov", _tower_DRCALO_NCerenkov, "tower_DRCALO_NCerenkov[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_iEta", _tower_DRCALO_iEta, "tower_DRCALO_iEta[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_iPhi", _tower_DRCALO_iPhi, "tower_DRCALO_iPhi[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_trueID", _tower_DRCALO_trueID, "tower_DRCALO_trueID[tower_DRCALO_N]/I");
  }
  if (_do_LFHCAL)
  {
    // towers LFHCALO
    _event_tree->Branch("tower_LFHCAL_N", &_nTowers_LFHCAL, "tower_LFHCAL_N/I");
    _event_tree->Branch("tower_LFHCAL_E", _tower_LFHCAL_E, "tower_LFHCAL_E[tower_LFHCAL_N]/F");
    _event_tree->Branch("tower_LFHCAL_iEta", _tower_LFHCAL_iEta, "tower_LFHCAL_iEta[tower_LFHCAL_N]/I");
    _event_tree->Branch("tower_LFHCAL_iPhi", _tower_LFHCAL_iPhi, "tower_LFHCAL_iPhi[tower_LFHCAL_N]/I");
    _event_tree->Branch("tower_LFHCAL_iL", _tower_LFHCAL_iL, "tower_LFHCAL_iL[tower_LFHCAL_N]/I");
    _event_tree->Branch("tower_LFHCAL_trueID", _tower_LFHCAL_trueID, "tower_LFHCAL_trueID[tower_LFHCAL_N]/I");
  }
  if (_do_FEMC)
  {
    // towers FEMC
    _event_tree->Branch("tower_FEMC_N", &_nTowers_FEMC, "tower_FEMC_N/I");
    _event_tree->Branch("tower_FEMC_E", _tower_FEMC_E, "tower_FEMC_E[tower_FEMC_N]/F");
    _event_tree->Branch("tower_FEMC_iEta", _tower_FEMC_iEta, "tower_FEMC_iEta[tower_FEMC_N]/I");
    _event_tree->Branch("tower_FEMC_iPhi", _tower_FEMC_iPhi, "tower_FEMC_iPhi[tower_FEMC_N]/I");
    _event_tree->Branch("tower_FEMC_trueID", _tower_FEMC_trueID, "tower_FEMC_trueID[tower_FEMC_N]/I");
    if (_do_CLUSTERS)
    {
      // clusters FEMC
      _event_tree->Branch("cluster_FEMC_N", &_nclusters_FEMC, "cluster_FEMC_N/I");
      _event_tree->Branch("cluster_FEMC_E", _cluster_FEMC_E, "cluster_FEMC_E[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_Eta", _cluster_FEMC_Eta, "cluster_FEMC_Eta[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_Phi", _cluster_FEMC_Phi, "cluster_FEMC_Phi[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_NTower", _cluster_FEMC_NTower, "cluster_FEMC_NTower[cluster_FEMC_N]/I");
      _event_tree->Branch("cluster_FEMC_trueID", _cluster_FEMC_trueID, "cluster_FEMC_trueID[cluster_FEMC_N]/I");
    }
  }
  if (_do_EEMC)
  {
    // towers EEMC
    _event_tree->Branch("tower_EEMC_N", &_nTowers_EEMC, "tower_EEMC_N/I");
    _event_tree->Branch("tower_EEMC_E", _tower_EEMC_E, "tower_EEMC_E[tower_EEMC_N]/F");
    _event_tree->Branch("tower_EEMC_iEta", _tower_EEMC_iEta, "tower_EEMC_iEta[tower_EEMC_N]/I");
    _event_tree->Branch("tower_EEMC_iPhi", _tower_EEMC_iPhi, "tower_EEMC_iPhi[tower_EEMC_N]/I");
    _event_tree->Branch("tower_EEMC_trueID", _tower_EEMC_trueID, "tower_EEMC_trueID[tower_EEMC_N]/I");
    if (_do_CLUSTERS)
    {
      // clusters EEMC
      _event_tree->Branch("cluster_EEMC_N", &_nclusters_EEMC, "cluster_EEMC_N/I");
      _event_tree->Branch("cluster_EEMC_E", _cluster_EEMC_E, "cluster_EEMC_E[cluster_EEMC_N]/F");
      _event_tree->Branch("cluster_EEMC_Eta", _cluster_EEMC_Eta, "cluster_EEMC_Eta[cluster_EEMC_N]/F");
      _event_tree->Branch("cluster_EEMC_Phi", _cluster_EEMC_Phi, "cluster_EEMC_Phi[cluster_EEMC_N]/F");
      _event_tree->Branch("cluster_EEMC_NTower", _cluster_EEMC_NTower, "cluster_EEMC_NTower[cluster_EEMC_N]/I");
      _event_tree->Branch("cluster_EEMC_trueID", _cluster_EEMC_trueID, "cluster_EEMC_trueID[cluster_EEMC_N]/I");
    }
  }
  
  if (_do_VERTEX)
  {
    // vertex
    _event_tree->Branch("vertex_x", &_vertex_x, "vertex_x/F");
    _event_tree->Branch("vertex_y", &_vertex_y, "vertex_y/F");
    _event_tree->Branch("vertex_z", &_vertex_z, "vertex_z/F");
    _event_tree->Branch("vertex_NCont", &_vertex_NCont, "vertex_NCont/I");
    _event_tree->Branch("vertex_true_x", &_vertex_true_x, "vertex_true_x/F");
    _event_tree->Branch("vertex_true_y", &_vertex_true_y, "vertex_true_y/F");
    _event_tree->Branch("vertex_true_z", &_vertex_true_z, "vertex_true_z/F");
  }
  if (_do_MCPARTICLES)
  {
    // MC particles
    _event_tree->Branch("nMCPart", &_nMCPart, "nMCPart/I");
    _event_tree->Branch("mcpart_ID", _mcpart_ID, "mcpart_ID[nMCPart]/I");
    _event_tree->Branch("mcpart_ID_parent", _mcpart_ID_parent, "mcpart_ID_parent[nMCPart]/I");
    _event_tree->Branch("mcpart_PDG", _mcpart_PDG, "mcpart_PDG[nMCPart]/I");
    _event_tree->Branch("mcpart_E", _mcpart_E, "mcpart_E[nMCPart]/F");
    _event_tree->Branch("mcpart_px", _mcpart_px, "mcpart_px[nMCPart]/F");
    _event_tree->Branch("mcpart_py", _mcpart_py, "mcpart_py[nMCPart]/F");
    _event_tree->Branch("mcpart_pz", _mcpart_pz, "mcpart_pz[nMCPart]/F");
    _event_tree->Branch("mcpart_x", _mcpart_x, "mcpart_x[nMCPart]/F");
    _event_tree->Branch("mcpart_y", _mcpart_y, "mcpart_y[nMCPart]/F");
    _event_tree->Branch("mcpart_z", _mcpart_z, "mcpart_z[nMCPart]/F");
    _event_tree->Branch("mcpart_BCID", _mcpart_BCID, "mcpart_BCID[nMCPart]/I");
  }
  if (_do_HEPMC)
  {
    // MC particles
    _event_tree->Branch("nHepmcp", &_nHepmcp, "nHepmcp/I");
    _event_tree->Branch("hepmcp_procid", &_hepmcp_procid, "hepmcp_procid/I");
    _event_tree->Branch("hepmcp_x1", &_hepmcp_x1, "hepmcp_x1/F");
    _event_tree->Branch("hepmcp_x2", &_hepmcp_x2, "hepmcp_x2/F");
    _event_tree->Branch("hepmcp_Q2", &_hepmcp_Q2, "hepmcp_Q2/F");
    _event_tree->Branch("hepmcp_vtx_x", &_hepmcp_vtx_x, "hepmcp_vtx_x/F");
    _event_tree->Branch("hepmcp_vtx_y", &_hepmcp_vtx_y, "hepmcp_vtx_y/F");
    _event_tree->Branch("hepmcp_vtx_z", &_hepmcp_vtx_z, "hepmcp_vtx_z/F");
    _event_tree->Branch("hepmcp_vtx_t", &_hepmcp_vtx_t, "hepmcp_vtx_t/F");

    //    _event_tree->Branch("hepmcp_ID_parent", _hepmcp_ID_parent, "hepmcp_ID_parent[nHepmcp]/F");
    _event_tree->Branch("hepmcp_status", _hepmcp_status, "hepmcp_status[nHepmcp]/I");
    _event_tree->Branch("hepmcp_PDG", _hepmcp_PDG, "hepmcp_PDG[nHepmcp]/I");
    _event_tree->Branch("hepmcp_E", _hepmcp_E, "hepmcp_E[nHepmcp]/F");
    _event_tree->Branch("hepmcp_px", _hepmcp_px, "hepmcp_px[nHepmcp]/F");
    _event_tree->Branch("hepmcp_py", _hepmcp_py, "hepmcp_py[nHepmcp]/F");
    _event_tree->Branch("hepmcp_pz", _hepmcp_pz, "hepmcp_pz[nHepmcp]/F");
    _event_tree->Branch("hepmcp_BCID", _hepmcp_BCID, "hepmcp_BCID[nHepmcp]/I");
    _event_tree->Branch("hepmcp_m1", _hepmcp_m1, "hepmcp_m1[nHepmcp]/I");
    _event_tree->Branch("hepmcp_m2", _hepmcp_m2, "hepmcp_m2[nHepmcp]/I");
  }

  if (_do_GEOMETRY)
  {
    _tfile_geometry = new TFile("geometry.root", "RECREATE");

    _geometry_tree = new TTree("geometry_tree", "geometry_tree");
    // tracks and hits
    _geometry_tree->Branch("calo", &_calo_ID, "nHits/I");
    _geometry_tree->Branch("calo_towers_N", &_calo_towers_N, "calo_towers_N/I");
    _geometry_tree->Branch("calo_towers_iEta", _calo_towers_iEta, "calo_towers_iEta[calo_towers_N]/I");
    _geometry_tree->Branch("calo_towers_iPhi", _calo_towers_iPhi, "calo_towers_iPhi[calo_towers_N]/I");
    _geometry_tree->Branch("calo_towers_iL", _calo_towers_iL, "calo_towers_iL[calo_towers_N]/I");
    _geometry_tree->Branch("calo_towers_Eta", _calo_towers_Eta, "calo_towers_Eta[calo_towers_N]/F");
    _geometry_tree->Branch("calo_towers_Phi", _calo_towers_Phi, "calo_towers_Phi[calo_towers_N]/F");
    _geometry_tree->Branch("calo_towers_x", _calo_towers_x, "calo_towers_x[calo_towers_N]/F");
    _geometry_tree->Branch("calo_towers_y", _calo_towers_y, "calo_towers_y[calo_towers_N]/F");
    _geometry_tree->Branch("calo_towers_z", _calo_towers_z, "calo_towers_z[calo_towers_N]/F");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventEvaluatorEIC::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
    std::cout << "entered process_event" << endl;

  if (_do_BECAL)
  {
    if (!_caloevalstackBECAL)
    {
      _caloevalstackBECAL = new CaloEvalStack(topNode, "BECAL");
      _caloevalstackBECAL->set_strict(_strict);
      _caloevalstackBECAL->set_verbosity(Verbosity() + 1);
    }
    else
      _caloevalstackBECAL->next_event(topNode);
  }
  if (_do_HCALIN)
  {
    if (!_caloevalstackHCALIN)
    {
      _caloevalstackHCALIN = new CaloEvalStack(topNode, "HCALIN");
      _caloevalstackHCALIN->set_strict(_strict);
      _caloevalstackHCALIN->set_verbosity(Verbosity() + 1);
    }
    else
      _caloevalstackHCALIN->next_event(topNode);
  }
  if (_do_HCALOUT)
  {
    if (!_caloevalstackHCALOUT)
    {
      _caloevalstackHCALOUT = new CaloEvalStack(topNode, "HCALOUT");
      _caloevalstackHCALOUT->set_strict(_strict);
      _caloevalstackHCALOUT->set_verbosity(Verbosity() + 1);
    }
    else
      _caloevalstackHCALOUT->next_event(topNode);
  }
  if (_do_EHCAL)
  {
    if (!_caloevalstackEHCAL)
    {
      _caloevalstackEHCAL = new CaloEvalStack(topNode, "EHCAL");
      _caloevalstackEHCAL->set_strict(_strict);
      _caloevalstackEHCAL->set_verbosity(Verbosity() + 1);
    }
    else
      _caloevalstackEHCAL->next_event(topNode);
  }
  if (_do_DRCALO)
  {
    if (!_caloevalstackDRCALO)
    {
      _caloevalstackDRCALO = new CaloEvalStack(topNode, "DRCALO");
      _caloevalstackDRCALO->set_strict(_strict);
      _caloevalstackDRCALO->set_verbosity(Verbosity() + 1);
    }
    else
      _caloevalstackDRCALO->next_event(topNode);
  }
  if (_do_LFHCAL)
  {
    if (!_caloevalstackLFHCAL)
    {
      _caloevalstackLFHCAL = new CaloEvalStack(topNode, "LFHCAL");
      _caloevalstackLFHCAL->set_strict(_strict);
      _caloevalstackLFHCAL->set_verbosity(Verbosity() + 1);
    }
    else
      _caloevalstackLFHCAL->next_event(topNode);
  }

  if (_do_FEMC)
  {
    if (!_caloevalstackFEMC)
    {
      _caloevalstackFEMC = new CaloEvalStack(topNode, "FEMC");
      _caloevalstackFEMC->set_strict(_strict);
      _caloevalstackFEMC->set_verbosity(Verbosity() + 1);
    }
    else
      _caloevalstackFEMC->next_event(topNode);
  }
  if (_do_EEMC)
  {
    if (!_caloevalstackEEMC)
    {
      _caloevalstackEEMC = new CaloEvalStack(topNode, "EEMC");
      _caloevalstackEEMC->set_strict(_strict);
      _caloevalstackEEMC->set_verbosity(Verbosity() + 1);
    }
    else
      _caloevalstackEEMC->next_event(topNode);
  }

  if (Verbosity() > 0)
    std::cout << "loaded evalstack" << endl;

  // fill the Evaluator Tree
  fillOutputNtuples(topNode);
  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

void EventEvaluatorEIC::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 2)
    std::cout << "EventEvaluatorEIC::fillOutputNtuples() entered" << endl;

  //----------------------
  // fill the Event Tree
  //----------------------

  //----------------------
  // Event level info
  //---------------------
  // Extract weight info from the stored HepMC event.
  if (_do_store_event_info)
  {
    PHHepMCGenEventMap* hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    if (hepmceventmap)
    {
      if (Verbosity() > 0)
        std::cout << "saving event level info" << endl;

      for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
           eventIter != hepmceventmap->end();
           ++eventIter)
      {
        PHHepMCGenEvent* hepmcevent = eventIter->second;

        if (hepmcevent)
        {
          HepMC::GenEvent* truthevent = hepmcevent->getEvent();
          if (!truthevent)
          {
            std::cout << PHWHERE
                 << "no evt pointer under phhepmvgeneventmap found "
                 << endl;
            return;
          }

          auto xsec = truthevent->cross_section();
          if (xsec)
            _cross_section = xsec->cross_section();

          // Only fill the event weight if available.
          // The overall event weight will be stored in the last entry in the vector.
          auto weights = truthevent->weights();
          if (weights.size() > 0)
            _event_weight = weights[weights.size() - 1];
        }
      }
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " PHHepMCGenEventMap node (for event level info) not found on node tree" << endl;
      return;
    }

    // Retrieve the number of generator accepted events
    // Following how this was implemented in PHPythia8
    PHNodeIterator iter(topNode);
    PHCompositeNode* sumNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
    if (!sumNode)
    {
      std::cout << PHWHERE << "RUN Node missing doing nothing" << endl;
      return;
    }
    auto* integralNode = findNode::getClass<PHGenIntegral>(sumNode, "PHGenIntegral");
    if (integralNode)
    {
      _n_generator_accepted = integralNode->get_N_Generator_Accepted_Event();
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " PHGenIntegral node (for n generator accepted) not found on node tree. Continuing" << endl;
    }
  }
  //----------------------
  //    VERTEX
  //----------------------
  SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (_do_VERTEX)
  {
    if (vertexmap)
    {
      if (!vertexmap->empty())
      {
        if (Verbosity() > 0)
          std::cout << "saving vertex" << endl;

        SvtxVertex* vertex = (vertexmap->begin())->second;
        _vertex_x = vertex->get_x();
        _vertex_y = vertex->get_y();
        _vertex_z = vertex->get_z();
        _vertex_NCont = vertex->size_tracks();
      }
      else
      {
        _vertex_x = 0.;
        _vertex_y = 0.;
        _vertex_z = 0.;
        _vertex_NCont = -1;
      }
    }
  }
  
  Int_t hitContainerIDs[200]    = {0};

  for (int iIndex = 0; iIndex < 200; ++iIndex){
    // you need to add your layer name here to be saved! This has to be done
    // as we do not want to save thousands of calorimeter hits!
    if (
          (GetProjectionNameFromIndex(iIndex).find("LFHCAL") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("HCALIN") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("HCALOUT") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("BECAL") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("FEMC") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("EHCAL") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("EEMCH") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("EEMC") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("TTL") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("LBLVTX") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("BARR") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("BST") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("FST") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("SVTX") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("EST") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("ZDCsurrogate") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("rpTruth") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("rpTruth2") != std::string::npos) ||  // needed for IP8
          (GetProjectionNameFromIndex(iIndex).find("offMomTruth") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("b0Truth") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("EGEM") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("FGEM") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("RWELL") != std::string::npos) 
      ){
      string nodename = "G4HIT_" + GetProjectionNameFromIndex(iIndex);
      PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (hits){
        hitContainerIDs[iIndex] = hits->GetID();
        if (Verbosity() > 3) std::cout << nodename << "\t" << hits->GetID() << "\t" << hitContainerIDs[iIndex]  << endl;
      }
    }
  }
  
  //----------------------
  //    HITS
  //----------------------
  if (_do_HITS)
  {
    if (Verbosity() > 0)
      std::cout << "saving hits" << endl;

    _nHitsLayers = 0;
    PHG4TruthInfoContainer* truthinfocontainerHits = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    for (int iIndex = 0; iIndex < 200; ++iIndex)
    {
      // you need to add your layer name here to be saved! This has to be done
      // as we do not want to save thousands of calorimeter hits!
      if (
          (_do_HITS_CALO && ( (GetProjectionNameFromIndex(iIndex).find("LFHCAL") != std::string::npos) ||
                              (GetProjectionNameFromIndex(iIndex).find("HCALIN") != std::string::npos) ||
                              (GetProjectionNameFromIndex(iIndex).find("HCALOUT") != std::string::npos) ||
                              (GetProjectionNameFromIndex(iIndex).find("BECAL") != std::string::npos) ||
                              (GetProjectionNameFromIndex(iIndex).find("FEMC") != std::string::npos) ||
                              (GetProjectionNameFromIndex(iIndex).find("EHCAL") != std::string::npos) ||
                              (GetProjectionNameFromIndex(iIndex).find("EEMCH") != std::string::npos) ||
                              (GetProjectionNameFromIndex(iIndex).find("EEMC") != std::string::npos))) ||
          (GetProjectionNameFromIndex(iIndex).find("TTL") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("BST") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("LBLVTX") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("BARR") != std::string::npos) ||
          // (GetProjectionNameFromIndex(iIndex).find("TrackingService") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("FST") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("SVTX") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("EST") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("ZDCsurrogate") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("rpTruth") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("rpTruth2") != std::string::npos) ||  // needed for IP8
          (GetProjectionNameFromIndex(iIndex).find("offMomTruth") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("b0Truth") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("EGEM") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("FGEM") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("RWELL") != std::string::npos) ||
          (((GetProjectionNameFromIndex(iIndex).find("BH_1") != std::string::npos) || (GetProjectionNameFromIndex(iIndex).find("BH_FORWARD_PLUS") != std::string::npos) || (GetProjectionNameFromIndex(iIndex).find("BH_FORWARD_NEG") != std::string::npos)) && _do_BLACKHOLE))
      {
        string nodename = "G4HIT_" + GetProjectionNameFromIndex(iIndex);
        PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
        if (hits)
        {
          if (Verbosity() > 1)
            std::cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;
       
          PHG4HitContainer::ConstRange hit_range = hits->getHits();
          for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
          {
            if (Verbosity() > 1)
              std::cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;
            if (_nHitsLayers > _maxNHits)
            {
              std::cout << __PRETTY_FUNCTION__ << " exceededed maximum hit array size! Please check where these hits come from!" << endl;
              break;
            }

            _hits_x[_nHitsLayers] = hit_iter->second->get_x(0);
            _hits_y[_nHitsLayers] = hit_iter->second->get_y(0);
            _hits_z[_nHitsLayers] = hit_iter->second->get_z(0);
            _hits_x2[_nHitsLayers] = hit_iter->second->get_x(1);
            _hits_y2[_nHitsLayers] = hit_iter->second->get_y(1);
            _hits_z2[_nHitsLayers] = hit_iter->second->get_z(1);
            _hits_t[_nHitsLayers] = hit_iter->second->get_t(0);
            _hits_edep[_nHitsLayers] = hit_iter->second->get_edep();
            _hits_lightyield[_nHitsLayers] = hit_iter->second->get_light_yield();
            _hits_isAbsorber[_nHitsLayers] = 0;
            _hits_layerID[_nHitsLayers] = iIndex;
            // TODO possibly add layer information to output (e.g. BST layers)
            // _hits_detectorLayer[_nHitsLayers] = hit_iter->second->get_index_i();
              // std::cout << "i " << hit_iter->second->get_index_i() << std::endl;
            // std::cout << "i " << hit_iter->second->get_index_i() << "\tj " <<hit_iter->second->get_index_j() << "\tk " <<hit_iter->second->get_index_k() << "\tl " << hit_iter->second->get_index_l() << "\tsens_x "<< hit_iter->second->get_strip_z_index()<< "\tsens_y "<< hit_iter->second->get_strip_y_index()   << endl;
            if (truthinfocontainerHits)
            {
              PHG4Particle* particle = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());

              if (particle->get_parent_id() != 0)
              {
                PHG4Particle* g4particleMother = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());
                int mcSteps = 0;
                while (g4particleMother->get_parent_id() != 0)
                {
                  g4particleMother = truthinfocontainerHits->GetParticle(g4particleMother->get_parent_id());
                  if (g4particleMother == NULL) break;
                  mcSteps += 1;
                }
                if (mcSteps <= _depth_MCstack)
                  _hits_trueID[_nHitsLayers] = hit_iter->second->get_trkid();
                else
                {
                  PHG4Particle* g4particleMother2 = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());
                  int mcSteps2 = 0;
                  while (g4particleMother2->get_parent_id() != 0 && (mcSteps2 < (mcSteps - _depth_MCstack + 1)))
                  {
                    g4particleMother2 = truthinfocontainerHits->GetParticle(g4particleMother2->get_parent_id());
                    if (g4particleMother2 == NULL)
                      break;
                    else
                    {
                      _hits_trueID[_nHitsLayers] = g4particleMother2->get_parent_id();
                      mcSteps2 += 1;
                    }
                  }
                }
              }
              else
                _hits_trueID[_nHitsLayers] = hit_iter->second->get_trkid();
            }
            _nHitsLayers++;
          }
          if (Verbosity() > 0)
            std::cout << "saved\t" << _nHitsLayers << "\thits for " << GetProjectionNameFromIndex(iIndex) << endl;
        }
        else
        {
          if (Verbosity() > 0)
            std::cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;
          continue;
        }
        if(_do_HITS_ABSORBER){
          string nodenameAbs = "G4HIT_ABSORBER_" + GetProjectionNameFromIndex(iIndex);
          PHG4HitContainer* hitsAbs = findNode::getClass<PHG4HitContainer>(topNode, nodenameAbs);
          if (hitsAbs)
          {
            if (Verbosity() > 1)
              std::cout << __PRETTY_FUNCTION__ << " absorber number of hits: " << hitsAbs->size() << endl;

            PHG4HitContainer::ConstRange hit_range = hitsAbs->getHits();
            for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
            {
              if (Verbosity() > 1)
                std::cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;
              if (_nHitsLayers > _maxNHits)
              {
                std::cout << __PRETTY_FUNCTION__ << " exceededed maximum hit array size! Please check where these hits come from!" << endl;
                break;
              }

              _hits_x[_nHitsLayers] = hit_iter->second->get_x(0);
              _hits_y[_nHitsLayers] = hit_iter->second->get_y(0);
              _hits_z[_nHitsLayers] = hit_iter->second->get_z(0);
              _hits_x2[_nHitsLayers] = hit_iter->second->get_x(1);
              _hits_y2[_nHitsLayers] = hit_iter->second->get_y(1);
              _hits_z2[_nHitsLayers] = hit_iter->second->get_z(1);
              _hits_t[_nHitsLayers] = hit_iter->second->get_t(0);
              _hits_edep[_nHitsLayers] = hit_iter->second->get_edep();
              _hits_lightyield[_nHitsLayers] = hit_iter->second->get_light_yield();
              _hits_isAbsorber[_nHitsLayers] = 1;
              _hits_layerID[_nHitsLayers] = iIndex;
//               std::cout << "i " << hit_iter->second->get_index_i() << "\tj " <<hit_iter->second->get_index_j() << "\tk " <<hit_iter->second->get_index_k() << "\tl " << hit_iter->second->get_index_l() << "\tsens_x "<< hit_iter->second->get_strip_z_index()<< "\tsens_y "<< hit_iter->second->get_strip_y_index()   << endl;
              if (truthinfocontainerHits)
              {
                PHG4Particle* particle = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());

                if (particle->get_parent_id() != 0)
                {
                  PHG4Particle* g4particleMother = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());
                  int mcSteps = 0;
                  while (g4particleMother->get_parent_id() != 0)
                  {
                    g4particleMother = truthinfocontainerHits->GetParticle(g4particleMother->get_parent_id());
                    if (g4particleMother == NULL) break;
                    mcSteps += 1;
                  }
                  if (mcSteps <= _depth_MCstack)
                    _hits_trueID[_nHitsLayers] = hit_iter->second->get_trkid();
                  else
                  {
                    PHG4Particle* g4particleMother2 = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());
                    int mcSteps2 = 0;
                    while (g4particleMother2->get_parent_id() != 0 && (mcSteps2 < (mcSteps - _depth_MCstack + 1)))
                    {
                      g4particleMother2 = truthinfocontainerHits->GetParticle(g4particleMother2->get_parent_id());
                      if (g4particleMother2 == NULL)
                        break;
                      else
                      {
                        _hits_trueID[_nHitsLayers] = g4particleMother2->get_parent_id();
                        mcSteps2 += 1;
                      }
                    }
                  }
                }
                else
                  _hits_trueID[_nHitsLayers] = hit_iter->second->get_trkid();
              }
              _nHitsLayers++;
            }
            if (Verbosity() > 0)
              std::cout << "saved\t" << _nHitsLayers << "\thits for absorber of " << GetProjectionNameFromIndex(iIndex) << endl;
          }
          else
          {
            if (Verbosity() > 0)
              std::cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;
            continue;
          }
        }
      }
    }
  }
  //----------------------
  //    TOWERS BECAL
  //----------------------
  if (_do_BECAL)
  {
    CaloRawTowerEval* towerevalBECAL = _caloevalstackBECAL->get_rawtower_eval();
    _nTowers_BECAL = 0;
    string towernodeBECAL = "TOWER_CALIB_BECAL";
    RawTowerContainer* towersBECAL = findNode::getClass<RawTowerContainer>(topNode, towernodeBECAL.c_str());

    if (Verbosity() > 0)
    {
      RawTowerContainer* towersBECAL1 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_RAW_BECAL");
      RawTowerContainer* towersBECAL2 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_SIM_BECAL");
      std::cout << "BECAL sim: " << towersBECAL2->size() << endl;
      std::cout << "BECAL raw: " << towersBECAL1->size() << endl;
      std::cout << "BECAL calib: " << towersBECAL->size() << endl;
    }
    if (towersBECAL)
    {
      if (Verbosity() > 0)
        std::cout << "saving BECAL towers" << endl;

      string towergeomnodeBECAL = "TOWERGEOM_BECAL";
      RawTowerGeomContainer* towergeomBECAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeBECAL.c_str());
      if (towergeomBECAL)
      {
        if (_do_GEOMETRY && !_geometry_done[kBECAL])
        {
          // std::ostream *fout= new ofstream("test.list");
          RawTowerGeomContainer::ConstRange all_towers = towergeomBECAL->get_tower_geometries();
          // *fout << "Calorimeter ID: " << towergeomBECAL->get_calorimeter_id() << endl;
          // *fout << "size: " << towergeomBECAL->size() << endl;
          // towergeomBECAL->identify(*fout);
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
               it != all_towers.second; ++it)
          {
            _calo_ID = kBECAL;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_iL[_calo_towers_N] = -1;
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            // it->second->identify(*fout);
            _calo_towers_N++;
          }
          _geometry_done[kBECAL] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
          // _calo_ID = kBECAL;
          // _calo_ID = kBECAL;
        }

        RawTowerContainer::ConstRange begin_end = towersBECAL->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold[kBECAL]) continue;
            _tower_BECAL_iEta[_nTowers_BECAL] = tower->get_bineta();
            _tower_BECAL_iPhi[_nTowers_BECAL] = tower->get_binphi();
            _tower_BECAL_E[_nTowers_BECAL] = tower->get_energy();

            PHG4Particle* primary = towerevalBECAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_BECAL_trueID[_nTowers_BECAL] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalBECAL->get_energy_contribution(tower, primary);
            }
            else
              _tower_BECAL_trueID[_nTowers_BECAL] = -10;

            _nTowers_BECAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << " ERROR: Can't find " << towergeomnodeBECAL << endl;
        // return;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nTowers_BECAL << "\tBECAL towers" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " ERROR: Can't find " << towernodeBECAL << endl;
      // return;
    }
  }
  //----------------------
  //    TOWERS HCALIN
  //----------------------
  if (_do_HCALIN)
  {
    CaloRawTowerEval* towerevalHCALIN = _caloevalstackHCALIN->get_rawtower_eval();
    _nTowers_HCALIN = 0;
    string towernodeHCALIN = "TOWER_CALIB_HCALIN";
    RawTowerContainer* towersHCALIN = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALIN.c_str());
    if (towersHCALIN)
    {
      if (Verbosity() > 0)
        std::cout << "saving HCAL towers" << endl;

      string towergeomnodeHCALIN = "TOWERGEOM_HCALIN";
      RawTowerGeomContainer* towergeomHCALIN = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALIN.c_str());
      if (towergeomHCALIN)
      {
        if (_do_GEOMETRY && !_geometry_done[kHCALIN])
        {
          RawTowerGeomContainer::ConstRange all_towers = towergeomHCALIN->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
               it != all_towers.second; ++it)
          {
            _calo_ID = kHCALIN;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_iL[_calo_towers_N] = -1;
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kHCALIN] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
        RawTowerContainer::ConstRange begin_end = towersHCALIN->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold[kHCALIN]) continue;
            _tower_HCALIN_iEta[_nTowers_HCALIN] = tower->get_bineta();
            _tower_HCALIN_iPhi[_nTowers_HCALIN] = tower->get_binphi();
            _tower_HCALIN_E[_nTowers_HCALIN] = tower->get_energy();

            PHG4Particle* primary = towerevalHCALIN->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_HCALIN_trueID[_nTowers_HCALIN] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalHCALIN->get_energy_contribution(tower, primary);
            }
            else
              _tower_HCALIN_trueID[_nTowers_HCALIN] = -10;
            _nTowers_HCALIN++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALIN << endl;
        // return;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nTowers_HCALIN << "\tHCALIN towers" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " ERROR: Can't find " << towernodeHCALIN << endl;
      // return;
    }
  }
  //----------------------
  //    TOWERS HCALOUT
  //----------------------
  if (_do_HCALOUT)
  {
    CaloRawTowerEval* towerevalHCALOUT = _caloevalstackHCALOUT->get_rawtower_eval();
    _nTowers_HCALOUT = 0;
    string towernodeHCALOUT = "TOWER_CALIB_HCALOUT";
    RawTowerContainer* towersHCALOUT = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALOUT.c_str());
    if (towersHCALOUT)
    {
      if (Verbosity() > 0)
        std::cout << "saving HCAL towers" << endl;
      
      string towergeomnodeHCALOUT = "TOWERGEOM_HCALOUT";
      RawTowerGeomContainer* towergeomHCALOUT = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALOUT.c_str());
      if (towergeomHCALOUT)
      {
        if (_do_GEOMETRY && !_geometry_done[kHCALOUT])
        {
          RawTowerGeomContainer::ConstRange all_towers = towergeomHCALOUT->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
               it != all_towers.second; ++it)
          {
            _calo_ID = kHCALOUT;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_iL[_calo_towers_N] = -1;
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kHCALOUT] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
        RawTowerContainer::ConstRange begin_end = towersHCALOUT->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold[kHCALOUT]) continue;
            _tower_HCALOUT_iEta[_nTowers_HCALOUT] = tower->get_bineta();
            _tower_HCALOUT_iPhi[_nTowers_HCALOUT] = tower->get_binphi();
            _tower_HCALOUT_E[_nTowers_HCALOUT] = tower->get_energy();

            PHG4Particle* primary = towerevalHCALOUT->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_HCALOUT_trueID[_nTowers_HCALOUT] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalHCALOUT->get_energy_contribution(tower, primary);
            }
            else
              _tower_HCALOUT_trueID[_nTowers_HCALOUT] = -10;
            _nTowers_HCALOUT++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALOUT << endl;
        // return;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nTowers_HCALOUT << "\tHCALOUT towers" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " ERROR: Can't find " << towernodeHCALOUT << endl;
      // return;
    }
  }
  //----------------------
  //    TOWERS EHCAL
  //----------------------
  if (_do_EHCAL)
  {
    CaloRawTowerEval* towerevalEHCAL = _caloevalstackEHCAL->get_rawtower_eval();
    _nTowers_EHCAL = 0;
    string towernodeEHCAL = "TOWER_CALIB_EHCAL";
    RawTowerContainer* towersEHCAL = findNode::getClass<RawTowerContainer>(topNode, towernodeEHCAL.c_str());
    if (towersEHCAL)
    {
      if (Verbosity() > 0)
        std::cout << "saving HCAL towers" << endl;
      string towergeomnodeEHCAL = "TOWERGEOM_EHCAL";
      RawTowerGeomContainer* towergeomEHCAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeEHCAL.c_str());
      if (towergeomEHCAL)
      {
        if (_do_GEOMETRY && !_geometry_done[kEHCAL])
        {
          RawTowerGeomContainer::ConstRange all_towers = towergeomEHCAL->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
               it != all_towers.second; ++it)
          {
            _calo_ID = kEHCAL;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_iL[_calo_towers_N] = -1;
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kEHCAL] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
        RawTowerContainer::ConstRange begin_end = towersEHCAL->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold[kEHCAL]) continue;
            _tower_EHCAL_iEta[_nTowers_EHCAL] = tower->get_bineta();
            _tower_EHCAL_iPhi[_nTowers_EHCAL] = tower->get_binphi();
            _tower_EHCAL_E[_nTowers_EHCAL] = tower->get_energy();

            PHG4Particle* primary = towerevalEHCAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_EHCAL_trueID[_nTowers_EHCAL] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalEHCAL->get_energy_contribution(tower, primary);
            }
            else
              _tower_EHCAL_trueID[_nTowers_EHCAL] = -10;
            _nTowers_EHCAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << " ERROR: Can't find " << towergeomnodeEHCAL << endl;
        // return;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nTowers_EHCAL << "\tEHCAL towers" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " ERROR: Can't find " << towernodeEHCAL << endl;
      // return;
    }
  }
  //----------------------
  //    TOWERS DRCALO
  //----------------------
  if (_do_DRCALO)
  {
    CaloRawTowerEval* towerevalDRCALO = _caloevalstackDRCALO->get_rawtower_eval();
    _nTowers_DRCALO = 0;
    string towernodeDRCALO = "TOWER_CALIB_DRCALO";
    RawTowerContainer* towersDRCALO = findNode::getClass<RawTowerContainer>(topNode, towernodeDRCALO.c_str());
    if (towersDRCALO)
    {
      if (Verbosity() > 0)
        std::cout << "saving DRCALO towers" << endl;
      string towergeomnodeDRCALO = "TOWERGEOM_DRCALO";
      RawTowerGeomContainer* towergeomDRCALO = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeDRCALO.c_str());
      if (towergeomDRCALO)
      {
        if (_do_GEOMETRY && !_geometry_done[kDRCALO])
        {
          RawTowerGeomContainer::ConstRange all_towers = towergeomDRCALO->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
               it != all_towers.second; ++it)
          {
            _calo_ID = kDRCALO;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_iL[_calo_towers_N] = -1;
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kDRCALO] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
        if (Verbosity() > 0)
          std::cout << "found DRCALO geom" << endl;
        RawTowerContainer::ConstRange begin_end = towersDRCALO->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          // if (dynamic_cast<RawTowerv2 *>(tower) == nullptr)
          // {
          //   std::cout << __PRETTY_FUNCTION__ << " : Fatal Error! "
          //       << "Expect RawTowerv2, but found this tower:";
          // tower->identify();
          //   // exit(1);
          // }
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold[kDRCALO]) continue;

            _tower_DRCALO_iEta[_nTowers_DRCALO] = tower->get_bineta();
            _tower_DRCALO_iPhi[_nTowers_DRCALO] = tower->get_binphi();
            _tower_DRCALO_E[_nTowers_DRCALO] = tower->get_energy();
            // std::cout << __LINE__ << endl;
            _tower_DRCALO_NScint[_nTowers_DRCALO] = tower->get_scint_gammas();
            _tower_DRCALO_NCerenkov[_nTowers_DRCALO] = tower->get_cerenkov_gammas();
            // std::cout << "sci gammas: " << tower->get_scint_gammas() << "\tcerenk gammas: " << tower->get_cerenkov_gammas() << endl;

            PHG4Particle* primary = towerevalDRCALO->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_DRCALO_trueID[_nTowers_DRCALO] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalDRCALO->get_energy_contribution(tower, primary);
            }
            else
              _tower_DRCALO_trueID[_nTowers_DRCALO] = -10;
            _nTowers_DRCALO++;
          }
        }
        if (Verbosity() > 0)
          std::cout << "finished DRCALO twr loop" << endl;
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << " ERROR: Can't find " << towergeomnodeDRCALO << endl;
        // return;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nTowers_DRCALO << "\tDRCALO towers" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " ERROR: Can't find " << towernodeDRCALO << endl;
      // return;
    }
  }
  //----------------------
  //    TOWERS LFHCAL
  //----------------------
  if (_do_LFHCAL)
  {
    CaloRawTowerEval* towerevalLFHCAL = _caloevalstackLFHCAL->get_rawtower_eval();
    _nTowers_LFHCAL = 0;
    string towernodeLFHCAL = "TOWER_CALIB_LFHCAL";
    RawTowerContainer* towersLFHCAL = findNode::getClass<RawTowerContainer>(topNode, towernodeLFHCAL.c_str());
    if (Verbosity() > 1)
      std::cout << "reading towers: " << towersLFHCAL->size() << endl;
    if (towersLFHCAL)
    {
      if (Verbosity() > 0)
      {
        std::cout << "saving LFHCAL towers" << endl;
      }
      string towergeomnodeLFHCAL = "TOWERGEOM_LFHCAL";
      RawTowerGeomContainer* towergeomLFHCAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeLFHCAL.c_str());
      if (towergeomLFHCAL)
      {
        if (_do_GEOMETRY && !_geometry_done[kLFHCAL])
        {
          RawTowerGeomContainer::ConstRange all_towers = towergeomLFHCAL->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
               it != all_towers.second; ++it)
          {
            _calo_ID = kLFHCAL;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_iL[_calo_towers_N] = it->second->get_binl();
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kLFHCAL] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
        if (Verbosity() > 0)
          std::cout << "found LFHCAL geom" << endl;

        RawTowerContainer::ConstRange begin_end = towersLFHCAL->getTowers();
        RawTowerContainer::ConstIterator rtiter;

        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold[kLFHCAL]) continue;
            if (Verbosity() > 1) std::cout << "\n event eval: \t" << tower->get_energy() << "\t ieta: " << tower->get_bineta() << "\t iphi: " << tower->get_binphi() << "\t iZ: " << tower->get_binl() << endl;
            _tower_LFHCAL_iEta[_nTowers_LFHCAL] = tower->get_bineta();
            _tower_LFHCAL_iPhi[_nTowers_LFHCAL] = tower->get_binphi();
            _tower_LFHCAL_iL[_nTowers_LFHCAL] = tower->get_binl();
            _tower_LFHCAL_E[_nTowers_LFHCAL] = tower->get_energy();
            PHG4Particle* primary = towerevalLFHCAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_LFHCAL_trueID[_nTowers_LFHCAL] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalLFHCAL->get_energy_contribution(tower, primary);
            }
            else
              _tower_LFHCAL_trueID[_nTowers_LFHCAL] = -10;
            _nTowers_LFHCAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << " ERROR: Can't find " << towergeomnodeLFHCAL << endl;
        // return;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nTowers_LFHCAL << "\tLFHCAL towers" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " ERROR: Can't find " << towernodeLFHCAL << endl;
      // return;
    }
  }

  //----------------------
  //    TOWERS FEMC
  //----------------------
  if (_do_FEMC)
  {
    CaloRawTowerEval* towerevalFEMC = _caloevalstackFEMC->get_rawtower_eval();
    _nTowers_FEMC = 0;
    string towernodeFEMC = "TOWER_CALIB_FEMC";
    RawTowerContainer* towersFEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeFEMC.c_str());
    if (towersFEMC)
    {
      if (Verbosity() > 0)
        std::cout << "saving EMC towers" << endl;
      string towergeomnodeFEMC = "TOWERGEOM_FEMC";
      RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeFEMC.c_str());
      if (towergeom)
      {
        if (_do_GEOMETRY && !_geometry_done[kFEMC])
        {
          RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
               it != all_towers.second; ++it)
          {
            _calo_ID = kFEMC;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_iL[_calo_towers_N] = -1;
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kFEMC] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
        RawTowerContainer::ConstRange begin_end = towersFEMC->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold[kFEMC]) continue;

            _tower_FEMC_iEta[_nTowers_FEMC] = tower->get_bineta();
            _tower_FEMC_iPhi[_nTowers_FEMC] = tower->get_binphi();
            _tower_FEMC_E[_nTowers_FEMC] = tower->get_energy();

            PHG4Particle* primary = towerevalFEMC->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_FEMC_trueID[_nTowers_FEMC] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalFEMC->get_energy_contribution(tower, primary);
            }
            else
              _tower_FEMC_trueID[_nTowers_FEMC] = -10;

            // PHG4TruthInfoContainer* truthinfocontainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
            // RawTower::ShowerConstRange shower_range = tower->get_g4showers();
            // for (RawTower::ShowerConstIterator iter = shower_range.first;
            //     iter != shower_range.second;
            //     ++iter)
            // {
            //   PHG4Shower* shower = truthinfocontainer->GetShower(iter->first);
            //   PHG4Particle* particleParent = truthinfocontainer->GetParticle(shower->get_parent_particle_id());
            //     if (particleParent)
            //     {
            //       std::cout << "\t\tcurr shower parent id: " << shower->get_parent_particle_id() << " ("<< _tower_FEMC_trueID[_nTowers_FEMC] << ")\tPDG " << particleParent->get_pid() << "\tedep: " << shower->get_edep() << endl;
            //     } else {
            //       std::cout << "\t\tcurr shower parent id: " << shower->get_parent_particle_id() << endl;
            //     }
            //   for (PHG4Shower::ParticleIdIter jter = shower->begin_g4particle_id();jter != shower->end_g4particle_id(); ++jter)
            //   {
            //     int g4particle_id = *jter;
            //     PHG4Particle* particle = truthinfocontainer->GetParticle(g4particle_id);
            //     if (particle)
            //     {
            //       // if(particle->get_e()>0.1*shower->get_edep()){
            //         // std::cout << "\t\t\tparticle ID in shower: " << particle->get_track_id() << "\tPDG " << particle->get_pid() << "\tenergy: " << particle->get_e()<< endl;
            //       // }
            //     }
            //   }
            // }
            _nTowers_FEMC++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << " ERROR: Can't find " << towergeomnodeFEMC << endl;
        // return;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nTowers_FEMC << "\tFEMC towers" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " ERROR: Can't find " << towernodeFEMC << endl;
      // return;
    }
  }
  //----------------------
  //    TOWERS EEMC
  //----------------------
  if (_do_EEMC)
  {
    CaloRawTowerEval* towerevalEEMC = _caloevalstackEEMC->get_rawtower_eval();
    _nTowers_EEMC = 0;
    string towernodeEEMC = "TOWER_CALIB_EEMC";
    RawTowerContainer* towersEEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeEEMC.c_str());
    if (towersEEMC)
    {
      if (Verbosity() > 0)
        std::cout << "saving EMC towers" << endl;
      
      string towergeomnodeEEMC = "TOWERGEOM_EEMC";
      RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeEEMC.c_str());
      if (towergeom)
      {
        if (_do_GEOMETRY && !_geometry_done[kEEMC])
        {
          RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
               it != all_towers.second; ++it)
          {
            _calo_ID = kEEMC;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_iL[_calo_towers_N] = -1;
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kEEMC] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
        RawTowerContainer::ConstRange begin_end = towersEEMC->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold[kEEMC]) continue;

            _tower_EEMC_iEta[_nTowers_EEMC] = tower->get_bineta();
            _tower_EEMC_iPhi[_nTowers_EEMC] = tower->get_binphi();
            _tower_EEMC_E[_nTowers_EEMC] = tower->get_energy();

            PHG4Particle* primary = towerevalEEMC->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_EEMC_trueID[_nTowers_EEMC] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalEEMC->get_energy_contribution(tower, primary);
            }
            else
              _tower_EEMC_trueID[_nTowers_EEMC] = -10;
            _nTowers_EEMC++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << " ERROR: Can't find " << towergeomnodeEEMC << endl;
        // return;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nTowers_EEMC << "\tEEMC towers" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " ERROR: Can't find " << towernodeEEMC << endl;
      // return;
    }
  }

  //------------------------
  // CLUSTERS HCALIN
  //------------------------
  if (_do_HCALIN && _do_CLUSTERS)
  {
    CaloRawClusterEval* clusterevalHCALIN = _caloevalstackHCALIN->get_rawcluster_eval();
    _nclusters_HCALIN = 0;
    if (Verbosity() > 1)
      std::cout << "CaloEvaluator::filling gcluster ntuple..." << endl;

    string clusternodeHCALIN = "CLUSTER_HCALIN";
    RawClusterContainer* clustersHCALIN = findNode::getClass<RawClusterContainer>(topNode, clusternodeHCALIN.c_str());
    if (clustersHCALIN)
    {
      // for every cluster
      for (const auto& iterator : clustersHCALIN->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold[kHCALIN]) continue;

        _cluster_HCALIN_E[_nclusters_HCALIN] = cluster->get_energy();
        _cluster_HCALIN_NTower[_nclusters_HCALIN] = cluster->getNTowers();
        _cluster_HCALIN_Phi[_nclusters_HCALIN] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_HCALIN_Eta[_nclusters_HCALIN] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_HCALIN_Eta[_nclusters_HCALIN] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
        }
        else
          _cluster_HCALIN_Eta[_nclusters_HCALIN] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
        
        PHG4Particle* primary = clusterevalHCALIN->max_truth_primary_particle_by_energy(cluster);

        if (primary)
          _cluster_HCALIN_trueID[_nclusters_HCALIN] = primary->get_track_id();
        else
          _cluster_HCALIN_trueID[_nclusters_HCALIN] = -10;

        _nclusters_HCALIN++;
      }
    }
    else
      std::cerr << PHWHERE << " ERROR: Can't find " << clusternodeHCALIN << endl;
      // return;
    if (Verbosity() > 0)
    {
      std::cout << "saved\t" << _nclusters_HCALIN << "\tHCALIN clusters" << endl;
    }
  }
  //------------------------
  // CLUSTERS HCALOUT
  //------------------------
  if (_do_HCALOUT && _do_CLUSTERS)
  {
    CaloRawClusterEval* clusterevalHCALOUT = _caloevalstackHCALOUT->get_rawcluster_eval();
    _nclusters_HCALOUT = 0;
    if (Verbosity() > 1)
    {
      std::cout << "CaloEvaluator::filling gcluster ntuple..." << endl;
    }

    string clusternodeHCALOUT = "CLUSTER_HCALOUT";
    RawClusterContainer* clustersHCALOUT = findNode::getClass<RawClusterContainer>(topNode, clusternodeHCALOUT.c_str());
    if (clustersHCALOUT)
    {
      // for every cluster
      for (const auto& iterator : clustersHCALOUT->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold[kHCALOUT]) continue;

        _cluster_HCALOUT_E[_nclusters_HCALOUT] = cluster->get_energy();
        _cluster_HCALOUT_NTower[_nclusters_HCALOUT] = cluster->getNTowers();
        _cluster_HCALOUT_Phi[_nclusters_HCALOUT] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_HCALOUT_Eta[_nclusters_HCALOUT] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_HCALOUT_Eta[_nclusters_HCALOUT] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
        }
        else
          _cluster_HCALOUT_Eta[_nclusters_HCALOUT] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
     
        PHG4Particle* primary = clusterevalHCALOUT->max_truth_primary_particle_by_energy(cluster);

        if (primary)
          _cluster_HCALOUT_trueID[_nclusters_HCALOUT] = primary->get_track_id();
        else
          _cluster_HCALOUT_trueID[_nclusters_HCALOUT] = -10;
        _nclusters_HCALOUT++;
      }
    }
    else
    {
      std::cerr << PHWHERE << " ERROR: Can't find " << clusternodeHCALOUT << endl;
      // return;
    }
    if (Verbosity() > 0)
    {
      std::cout << "saved\t" << _nclusters_HCALOUT << "\tHCALOUT clusters" << endl;
    }
  }
  //------------------------
  // CLUSTERS EHCAL
  //------------------------
  if (_do_EHCAL && _do_CLUSTERS)
  {
    CaloRawClusterEval* clusterevalEHCAL = _caloevalstackEHCAL->get_rawcluster_eval();
    _nclusters_EHCAL = 0;
    if (Verbosity() > 1)
      std::cout << "CaloEvaluator::filling gcluster ntuple..." << endl;
    
    string clusternodeEHCAL = "CLUSTER_EHCAL";
    RawClusterContainer* clustersEHCAL = findNode::getClass<RawClusterContainer>(topNode, clusternodeEHCAL.c_str());
    if (clustersEHCAL)
    {
      // for every cluster
      for (const auto& iterator : clustersEHCAL->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold[kEHCAL]) continue;

        _cluster_EHCAL_E[_nclusters_EHCAL] = cluster->get_energy();
        _cluster_EHCAL_NTower[_nclusters_EHCAL] = cluster->getNTowers();
        _cluster_EHCAL_Phi[_nclusters_EHCAL] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_EHCAL_Eta[_nclusters_EHCAL] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_EHCAL_Eta[_nclusters_EHCAL] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
        }
        else
          _cluster_EHCAL_Eta[_nclusters_EHCAL] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
    
        PHG4Particle* primary = clusterevalEHCAL->max_truth_primary_particle_by_energy(cluster);

        if (primary)
          _cluster_EHCAL_trueID[_nclusters_EHCAL] = primary->get_track_id();
        else
          _cluster_EHCAL_trueID[_nclusters_EHCAL] = -10;
    
        _nclusters_EHCAL++;
      }
    }
    else
      std::cerr << PHWHERE << " ERROR: Can't find " << clusternodeEHCAL << endl;
      // return;
    if (Verbosity() > 0)
      std::cout << "saved\t" << _nclusters_EHCAL << "\tEHCAL clusters" << endl;    
  }

  //------------------------
  // CLUSTERS FEMC
  //------------------------
  if (_do_FEMC && _do_CLUSTERS)
  {
    CaloRawClusterEval* clusterevalFEMC = _caloevalstackFEMC->get_rawcluster_eval();
    _nclusters_FEMC = 0;
    if (Verbosity() > 1)
      std::cout << "CaloEvaluator::filling gcluster ntuple..." << endl;
    
    string clusternodeFEMC = "CLUSTER_FEMC";
    RawClusterContainer* clustersFEMC = findNode::getClass<RawClusterContainer>(topNode, clusternodeFEMC.c_str());
    if (clustersFEMC)
    {
      // for every cluster
      for (const auto& iterator : clustersFEMC->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold[kFEMC]) continue;

        _cluster_FEMC_E[_nclusters_FEMC] = cluster->get_energy();
        _cluster_FEMC_NTower[_nclusters_FEMC] = cluster->getNTowers();
        _cluster_FEMC_Phi[_nclusters_FEMC] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_FEMC_Eta[_nclusters_FEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_FEMC_Eta[_nclusters_FEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
        }
        else
          _cluster_FEMC_Eta[_nclusters_FEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
       
        PHG4Particle* primary = clusterevalFEMC->max_truth_primary_particle_by_energy(cluster);

        if (primary)
          _cluster_FEMC_trueID[_nclusters_FEMC] = primary->get_track_id();
        else
          _cluster_FEMC_trueID[_nclusters_FEMC] = -10;
       
        _nclusters_FEMC++;
      }
    }
    else
      std::cerr << PHWHERE << " ERROR: Can't find " << clusternodeFEMC << endl;
      // return;
    if (Verbosity() > 0)
      std::cout << "saved\t" << _nclusters_FEMC << "\tFEMC clusters" << endl;
  }
  //------------------------
  // CLUSTERS EEMC
  //------------------------
  if (_do_EEMC && _do_CLUSTERS)
  {
    CaloRawClusterEval* clusterevalEEMC = _caloevalstackEEMC->get_rawcluster_eval();
    _nclusters_EEMC = 0;
    if (Verbosity() > 1)
      std::cout << "CaloEvaluator::filling gcluster ntuple..." << endl;

    string clusternodeEEMC = "CLUSTER_EEMC";
    RawClusterContainer* clustersEEMC = findNode::getClass<RawClusterContainer>(topNode, clusternodeEEMC.c_str());
    if (clustersEEMC)
    {
      // for every cluster
      for (const auto& iterator : clustersEEMC->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold[kEEMC]) continue;

        _cluster_EEMC_E[_nclusters_EEMC] = cluster->get_energy();
        _cluster_EEMC_NTower[_nclusters_EEMC] = cluster->getNTowers();
        _cluster_EEMC_Phi[_nclusters_EEMC] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_EEMC_Eta[_nclusters_EEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_EEMC_Eta[_nclusters_EEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));
        }
        else
          _cluster_EEMC_Eta[_nclusters_EEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));

        PHG4Particle* primary = clusterevalEEMC->max_truth_primary_particle_by_energy(cluster);

        if (primary)
          _cluster_EEMC_trueID[_nclusters_EEMC] = primary->get_track_id();
        else
          _cluster_EEMC_trueID[_nclusters_EEMC] = -10;

        _nclusters_EEMC++;
      }
    }
    else
      std::cerr << PHWHERE << " ERROR: Can't find " << clusternodeEEMC << endl;

    if (Verbosity() > 0)
      std::cout << "saved\t" << _nclusters_EEMC << "\tEEMC clusters" << endl;
  }

  //------------------------
  // TRACKS
  //------------------------
  if (_do_TRACKS)
  {
    _nTracks = 0;
    _nProjections = 0;

    EICPIDParticleContainer* pidcontainer(nullptr);

    if (_do_PID_LogLikelihood)
    {
      pidcontainer = findNode::getClass<EICPIDParticleContainer>(topNode, "EICPIDParticleMap");
      if (pidcontainer == nullptr)
        std::cout << __PRETTY_FUNCTION__ << " Error: missing EICPIDParticleMap while _do_PID_LogLikelihood = " << _do_PID_LogLikelihood << endl;
    }

    // Loop over track maps, identifiy each source.
    // Although this configuration is fixed here, it doesn't require multiple sources.
    // It will only store them if they're available.
    std::vector<std::pair<std::string, TrackSource_t>> trackMapPairs = {
        {"TrackMap", TrackSource_t::all},
        {"InnerTrackMap", TrackSource_t::inner},
        {"SiliconTrackMap", TrackSource_t::silicon},
        {"TTLTrackMap", TrackSource_t::ttl},
        {"DefaultTrackMap", TrackSource_t::defecce}
    };
    bool foundAtLeastOneTrackSource = false;
    for (const auto& trackMapInfo : trackMapPairs)
    {
      if (_nTracks >= _maxNTracks) break;

      SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, trackMapInfo.first);
      if (trackmap)
      {
        foundAtLeastOneTrackSource = true;
        int nTracksInASource = 0;
        if (Verbosity() > 0)
          std::cout << "saving tracks for track map: " << trackMapInfo.first << endl;
        for (SvtxTrackMap::ConstIter track_itr = trackmap->begin(); track_itr != trackmap->end(); track_itr++)
        {
          if (_nTracks >= _maxNTracks) break;

          SvtxTrack_FastSim* track = dynamic_cast<SvtxTrack_FastSim*>(track_itr->second);
          if (track)
          {
            _track_ID[_nTracks] = track->get_id();
            _track_charge[_nTracks] = track->get_charge();
            _track_px[_nTracks] = track->get_px();
            _track_py[_nTracks] = track->get_py();
            _track_pz[_nTracks] = track->get_pz();
            _track_x[_nTracks] = track->get_x();
            _track_y[_nTracks] = track->get_y();
            _track_z[_nTracks] = track->get_z();
            _track_ndf[_nTracks] = track->get_ndf();
            _track_chi2[_nTracks] = track->get_chisq();
            // Ideally, would be dca3d_xy and dca3d_z, but these don't seem to be calculated properly in the
            // current (June 2021) simulations (they return NaN). So we take dca (seems to be ~ the 3d distance)
            // and dca_2d (seems to be ~ the distance in the transverse plane).
            // The names of the branches are based on the method names.
            _track_dca[_nTracks] = static_cast<float>(track->get_dca());
            _track_dca_2d[_nTracks] = static_cast<float>(track->get_dca2d());
            _track_trueID[_nTracks] = track->get_truth_track_id();
            _track_source[_nTracks] = static_cast<unsigned short>(trackMapInfo.second);
            
            short nHits = 0;
            UInt_t hitencoded = 0;
            for (const auto &g4hit_id_hitset : track->g4hit_ids())
            {
              int iIndex = 0;
              // hit_iter->second->get_index_i()
              while ( g4hit_id_hitset.first!= 0 && g4hit_id_hitset.first != hitContainerIDs[iIndex] && iIndex < 200){
                iIndex++;
              }
              if ( GetExponentFromProjectionIndex(iIndex) != -1){
                //NOTE if the detector layer is already encoded, then we don't need to add it again
                // catches rare case of multiple hits in a layer
                if(hitencoded & (UInt_t)pow(2,GetExponentFromProjectionIndex(iIndex))){
                  if (Verbosity() > 3) std::cout << g4hit_id_hitset.first << "\t" << pow(2,GetExponentFromProjectionIndex(iIndex)) << " already encoded in hitencoded: " << hitencoded << endl;
                  continue;
                }
                hitencoded += pow(2,GetExponentFromProjectionIndex(iIndex));
                nHits++;
              }
              if (Verbosity() > 3) std::cout << g4hit_id_hitset.first << "\t" << GetExponentFromProjectionIndex(iIndex) << "\t" << GetProjectionNameFromIndex(iIndex) << endl;
            }
            if (Verbosity() > 3) std::cout << "Total hits\t" << nHits << "\t" << hitencoded <<  endl;
            _track_nHits[_nTracks]        = nHits;
            _track_hitsEncoded[_nTracks]  = hitencoded;
            
            
            if (_do_PROJECTIONS)
            {
              // find projections
              for (SvtxTrack::ConstStateIter trkstates = track->begin_states(); trkstates != track->end_states(); ++trkstates)
              {
                if (Verbosity() > 1)
                  std::cout << __PRETTY_FUNCTION__ << " processing " << trkstates->second->get_name() << endl;

                string trackStateName = trkstates->second->get_name();
                if (Verbosity() > 1)
                  std::cout << __PRETTY_FUNCTION__ << " found " << trkstates->second->get_name() << endl;

                int trackStateIndex = GetProjectionIndex(trackStateName);
                if (trackStateIndex > -1)
                {
                  // save true projection info to given branch
                  _track_TLP_x[_nProjections] = trkstates->second->get_pos(0);
                  _track_TLP_y[_nProjections] = trkstates->second->get_pos(1);
                  _track_TLP_z[_nProjections] = trkstates->second->get_pos(2);
                  _track_TLP_t[_nProjections] = trkstates->first;
                  _track_TLP_px[_nProjections] = trkstates->second->get_px();
                  _track_TLP_py[_nProjections] = trkstates->second->get_py();
                  _track_TLP_pz[_nProjections] = trkstates->second->get_pz();
                  _track_ProjLayer[_nProjections] = trackStateIndex;
                  _track_ProjTrackID[_nProjections] = _nTracks;

                  string nodename = "G4HIT_" + trkstates->second->get_name();
                  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
                  if (hits)
                  {
                    if (Verbosity() > 1)
                      std::cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;

                    PHG4HitContainer::ConstRange hit_range = hits->getHits();
                    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
                    {
                      if (Verbosity() > 2)
                        std::cout << __PRETTY_FUNCTION__ << " checking hit id " << hit_iter->second->get_trkid() << " against " << track->get_truth_track_id() << endl;

                      if (hit_iter->second->get_trkid() - track->get_truth_track_id() == 0)
                      {
                        if (Verbosity() > 2)
                          std::cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;

                        // save reco projection info to given branch
                        _track_TLP_true_x[_nProjections] = hit_iter->second->get_x(0);
                        _track_TLP_true_y[_nProjections] = hit_iter->second->get_y(0);
                        _track_TLP_true_z[_nProjections] = hit_iter->second->get_z(0);
                        _track_TLP_true_t[_nProjections] = hit_iter->second->get_t(0);
                      }
                    }
                  }
                  else
                  {
                    if (Verbosity() > 1)
                      std::cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;

                    _track_TLP_true_x[_nProjections] = -10000;
                    _track_TLP_true_y[_nProjections] = -10000;
                    _track_TLP_true_z[_nProjections] = -10000;
                    _track_TLP_true_t[_nProjections] = -10000;
                  }
                  _nProjections++;
                }
              }
            }  //             if (_do_PROJECTIONS)

            if (_do_PID_LogLikelihood)
            {
              // perform PID matching
              if ((trackMapInfo.second == TrackSource_t::all || trackMapInfo.second == TrackSource_t::defecce) and pidcontainer != nullptr)
              {
                // only do so for the TrackSource_t::all
                const EICPIDParticle* pid_particle =
                    pidcontainer->findEICPIDParticle(track->get_id());

                if (pid_particle)
                {
                  if ((unsigned int) _nTracks >= _track_pion_LL.size())
                  {
                    std::cout << __PRETTY_FUNCTION__
                         << " logical error _nTracks = " << _nTracks
                         << " logical error _track_pion_LL.size() = " << _track_pion_LL.size()
                         << endl;
                    exit(1);
                  }

                  _track_pion_LL[_nTracks] = pid_particle->get_SumLogLikelyhood(EICPIDDefs::PionCandiate);
                  _track_kaon_LL[_nTracks] = pid_particle->get_SumLogLikelyhood(EICPIDDefs::KaonCandiate);
                  _track_proton_LL[_nTracks] = pid_particle->get_SumLogLikelyhood(EICPIDDefs::ProtonCandiate);
                }
              }
            }  // if (_do_PID_LogLikelihood )

            _nTracks++;
            nTracksInASource++;
          }
          else
          {
            if (Verbosity() > 0)  //Verbosity()
            {
              std::cout << "PHG4TrackFastSimEval::fill_track_tree - ignore track that is not a SvtxTrack_FastSim:";
              track_itr->second->identify();
            }
            continue;
          }
        }
        if (Verbosity() > 0)
          std::cout << "saved\t" << nTracksInASource << "\ttracks from track map " << trackMapInfo.first << ". Total saved tracks: " << _nTracks << endl;
      }
      else
      {
        if (Verbosity() > 0)
          std::cout << PHWHERE << "SvtxTrackMap node with name '" << trackMapInfo.first << "' not found on node tree" << endl;
      }
    }
    if (foundAtLeastOneTrackSource == false)
    {
      std::cout << PHWHERE << "Requested tracks, but found no sources on node tree. Returning" << endl;
      return;
    }
  }
  //------------------------
  // MC PARTICLES
  //------------------------
  _nMCPart = 0;
  if (_do_MCPARTICLES)
  {
    PHG4TruthInfoContainer* truthinfocontainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfocontainer)
    {
      if (Verbosity() > 0)
        std::cout << "saving MC particles" << endl;

      //GetParticleRange for all particles
      //GetPrimaryParticleRange for primary particles
      PHG4TruthInfoContainer::ConstRange range = truthinfocontainer->GetParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first; truth_itr != range.second; ++truth_itr)
      {
        PHG4Particle* g4particle = truth_itr->second;
        if (!g4particle) continue;
        int mcSteps = 0;
        PHG4Particle* g4particleMother = truth_itr->second;
        if (g4particle->get_parent_id() != 0)
        {
          while (g4particleMother->get_parent_id() != 0)
          {
            g4particleMother = truthinfocontainer->GetParticle(g4particleMother->get_parent_id());
            if (g4particleMother == NULL) break;
            mcSteps += 1;
          }
        }
        if(_depth_MCstack==0){
          if (mcSteps > _depth_MCstack) continue;
        } else {
          if (g4particle->get_parent_id()<0 && mcSteps > _depth_MCstack) continue;
        }
        // evaluating true primary vertex
        if (_do_VERTEX && _nMCPart == 0)
        {
          PHG4VtxPoint* vtx = truthinfocontainer->GetVtx(g4particle->get_vtx_id());
          if (vtx)
          {
            _vertex_true_x = vtx->get_x();
            _vertex_true_y = vtx->get_y();
            _vertex_true_z = vtx->get_z();
          }
        }

        // in case of all MC particles, make restrictions on the secondary selection
        // if(g4particle->get_track_id()<0 && g4particle->get_e()<0.5) continue;
        // primary (g4particle->get_parent_id() == 0) selection via:
        // if(gtrackID < 0) continue;

        //using the e threshold also for the truth particles gets rid of all the low energy secondary particles
        if (g4particle->get_e() < _reco_e_thresholdMC) continue;

        _mcpart_ID[_nMCPart] = g4particle->get_track_id();
        _mcpart_ID_parent[_nMCPart] = g4particle->get_parent_id();
        _mcpart_PDG[_nMCPart] = g4particle->get_pid();
        _mcpart_E[_nMCPart] = g4particle->get_e();
        _mcpart_px[_nMCPart] = g4particle->get_px();
        _mcpart_py[_nMCPart] = g4particle->get_py();
        _mcpart_pz[_nMCPart] = g4particle->get_pz();
        PHG4VtxPoint* vtxtmp = truthinfocontainer->GetVtx(g4particle->get_vtx_id());
        if (vtxtmp)
        {
          _mcpart_x[_nMCPart] = vtxtmp->get_x();
          _mcpart_y[_nMCPart] = vtxtmp->get_y();
          _mcpart_z[_nMCPart] = vtxtmp->get_z();
        }
        //BCID added for G4Particle --  HEPMC particle matching
        _mcpart_BCID[_nMCPart] = g4particle->get_barcode();
        // TVector3 projvec(_mcpart_px[0],_mcpart_py[0],_mcpart_pz[0]);
        // float projeta = projvec.Eta();
        _nMCPart++;
      }
      if (Verbosity() > 0)
        std::cout << "saved\t" << _nMCPart << "\tMC particles" << endl;
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;
      return;
    }
  }

  //------------------------
  // HEPMC
  //------------------------
  _nHepmcp = 0;
  if (_do_HEPMC)
  {
    PHHepMCGenEventMap* hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    if (hepmceventmap)
    {
      if (Verbosity() > 0)
        std::cout << "saving HepMC output" << endl;

      if (Verbosity() > 0)
        hepmceventmap->Print();

      for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
           eventIter != hepmceventmap->end();
           ++eventIter)
      {
        PHHepMCGenEvent* hepmcevent = eventIter->second;

        if (hepmcevent)
        {
          HepMC::GenEvent* truthevent = hepmcevent->getEvent();
          if (!truthevent)
          {
            std::cout << PHWHERE
                 << "no evt pointer under phhepmvgeneventmap found "
                 << endl;
            return;
          }

          _hepmcp_vtx_x = hepmcevent->get_collision_vertex().x();
          _hepmcp_vtx_y = hepmcevent->get_collision_vertex().y();
          _hepmcp_vtx_z = hepmcevent->get_collision_vertex().z();
          _hepmcp_vtx_t = hepmcevent->get_collision_vertex().t();

          HepMC::PdfInfo* pdfinfo = truthevent->pdf_info();

          //     m_partid1 = pdfinfo->id1();
          // m_partid2 = pdfinfo->id2();
          _hepmcp_x1 = pdfinfo->x1();
          _hepmcp_x2 = pdfinfo->x2();
          _hepmcp_Q2 = pdfinfo->scalePDF();

          // m_mpi = truthevent->mpi();

          _hepmcp_procid = truthevent->signal_process_id();

          if (Verbosity() > 2)
          {
            std::cout << " Iterating over an event" << endl;
          }
          for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
               iter != truthevent->particles_end();
               ++iter)
          {
            _hepmcp_E[_nHepmcp] = (*iter)->momentum().e();
            _hepmcp_PDG[_nHepmcp] = (*iter)->pdg_id();
            _hepmcp_px[_nHepmcp] = (*iter)->momentum().px();
            _hepmcp_py[_nHepmcp] = (*iter)->momentum().py();
            _hepmcp_pz[_nHepmcp] = (*iter)->momentum().pz();
            _hepmcp_status[_nHepmcp] = (*iter)->status();
            _hepmcp_BCID[_nHepmcp] = (*iter)->barcode();
            _hepmcp_m2[_nHepmcp] = 0;
            _hepmcp_m1[_nHepmcp] = 0;
            if ((*iter)->production_vertex())
            {
              for (HepMC::GenVertex::particle_iterator mother = (*iter)->production_vertex()->particles_begin(HepMC::parents);
                   mother != (*iter)->production_vertex()->particles_end(HepMC::parents);
                   ++mother)
              {
                _hepmcp_m2[_nHepmcp] = (*mother)->barcode();
                if (_hepmcp_m1[_nHepmcp] == 0)
                  _hepmcp_m1[_nHepmcp] = (*mother)->barcode();
              }
            }
            if (Verbosity() > 2) 
              std::cout << "nHepmcp " << _nHepmcp << "\tPDG " << _hepmcp_PDG[_nHepmcp] << "\tEnergy " << _hepmcp_E[_nHepmcp] << "\tbarcode " << _hepmcp_BCID[_nHepmcp] << "\tMother1 " << _hepmcp_m1[_nHepmcp] << "\tMother2 " << _hepmcp_m2[_nHepmcp] << endl;
            _nHepmcp++;
          }
        }
      }
    }
    else
    {
      if (Verbosity() > 0)
        std::cout << PHWHERE << " PHHepMCGenEventMap node not found on node tree" << endl;
      return;
    }
  }  //hepmc

  _event_tree->Fill();

  if (Verbosity() > 0)
    std::cout << "Resetting buffer ..." << endl;

  resetBuffer();
  if (Verbosity() > 0)
    std::cout << "EventEvaluatorEIC buffer reset" << endl;
  return;
}

int EventEvaluatorEIC::End(PHCompositeNode* topNode)
{
  _tfile->cd();

  _event_tree->Write();

  _tfile->Close();

  delete _tfile;

  if (_do_GEOMETRY)
  {
    _tfile_geometry->cd();

    _geometry_tree->Write();

    _tfile_geometry->Close();

    delete _tfile_geometry;
  }
  if (Verbosity() > 0)
  {
    std::cout << "========================= " << Name() << "::End() ============================" << endl;
    std::cout << " " << _ievent << " events of output written to: " << _filename << endl;
    std::cout << "===========================================================================" << endl;
  }

  if (_caloevalstackBECAL) delete _caloevalstackBECAL;
  if (_caloevalstackHCALIN) delete _caloevalstackHCALIN;
  if (_caloevalstackHCALOUT) delete _caloevalstackHCALOUT;
  if (_caloevalstackEHCAL) delete _caloevalstackEHCAL;
  if (_caloevalstackDRCALO) delete _caloevalstackDRCALO;
  if (_caloevalstackLFHCAL) delete _caloevalstackLFHCAL;
  if (_caloevalstackFEMC) delete _caloevalstackFEMC;
  if (_caloevalstackEEMC) delete _caloevalstackEEMC;
 
  return Fun4AllReturnCodes::EVENT_OK;
}

int EventEvaluatorEIC::GetProjectionIndex(std::string projname)
{
  if (projname.find("FTTL_0") != std::string::npos)
    return 0;
  else if (projname.find("FTTL_1") != std::string::npos)
    return 1;
  else if (projname.find("FTTL_2") != std::string::npos)
    return 2;
  else if (projname.find("ETTL_0") != std::string::npos)
    return 3;
  else if (projname.find("ETTL_1") != std::string::npos)
    return 4;
  else if (projname.find("LFHCAL") != std::string::npos)
    return 5;
  else if (projname.find("FEMC") != std::string::npos)
    return 6;
  else if (projname.find("CTTL_0") != std::string::npos)
    return 7;
  else if (projname.find("CTTL_1") != std::string::npos)
    return 8;

  else if (projname.find("hpDIRC") != std::string::npos)
    return 35;
  else if (projname.find("mRICH_Plane") != std::string::npos)
    return 36;
  else if (projname.find("dRICH_Plane") != std::string::npos)
    return 37;

  else if (projname.find("FST_0") != std::string::npos)
    return 50;
  else if (projname.find("FST_1") != std::string::npos)
    return 51;
  else if (projname.find("FST_2") != std::string::npos)
    return 52;
  else if (projname.find("FST_3") != std::string::npos)
    return 53;
  else if (projname.find("FST_4") != std::string::npos)
    return 54;
  else if (projname.find("FST_5") != std::string::npos)
    return 55;
  else if (projname.find("FST_6") != std::string::npos)
    return 56;

  else if (projname.find("EHCAL") != std::string::npos)
    return 60;
  else if (projname.find("EEMC") != std::string::npos)
    return 61;
  else if (projname.find("HCALIN") != std::string::npos)
    return 62;
  else if (projname.find("HCALOUT") != std::string::npos)
    return 63;
  else if (projname.find("BECAL") != std::string::npos)
    return 66;
 
  else if (projname.find("ZDCsurrogate") != std::string::npos)
    return 70;
  else if (projname.find("rpTruth") != std::string::npos)
    return 71;
  else if (projname.find("b0Truth") != std::string::npos)
    return 72;
  else if (projname.find("rpTruth2") != std::string::npos)
    return 73;
  else if (projname.find("offMomTruth") != std::string::npos)
    return 74;

  else if (projname.find("TrackingService") != std::string::npos)
    return 79;

  else if (projname.find("BH_1") != std::string::npos)
    return 90;
  else if (projname.find("BH_FORWARD_PLUS") != std::string::npos)
    return 91;
  else if (projname.find("BH_FORWARD_NEG") != std::string::npos)
    return 92;
 
  else if (projname.find("RWELL_0") != std::string::npos)
    return 110;
  else if (projname.find("RWELL_1") != std::string::npos)
    return 111;
  else if (projname.find("RWELL_2") != std::string::npos)
    return 112;


  else if (projname.find("EST_0") != std::string::npos)
    return 113;
  else if (projname.find("EST_1") != std::string::npos)
    return 114;
  else if (projname.find("EST_2") != std::string::npos)
    return 115;
  else if (projname.find("EST_3") != std::string::npos)
    return 116;
  else if (projname.find("EST_4") != std::string::npos)
    return 117;
  else if (projname.find("EST_5") != std::string::npos)
    return 118;
  else if (projname.find("EST_6") != std::string::npos)
    return 119;

  else if (projname.find("FGEM_0") != std::string::npos)
    return 120;
  else if (projname.find("FGEM_1") != std::string::npos)
    return 121;

  else if (projname.find("EGEM_0") != std::string::npos)
    return 130;
  else if (projname.find("EGEM_1") != std::string::npos)
    return 131;


  else if (projname.find("LFHCAL_0") != std::string::npos)
    return 140;
  else if (projname.find("LFHCAL_1") != std::string::npos)
    return 141;
  else if (projname.find("LFHCAL_2") != std::string::npos)
    return 142;
  else if (projname.find("LFHCAL_3") != std::string::npos)
    return 143;
  else if (projname.find("LFHCAL_4") != std::string::npos)
    return 144;
  else if (projname.find("LFHCAL_5") != std::string::npos)
    return 145;
  else if (projname.find("LFHCAL_6") != std::string::npos)
    return 146;
  else if (projname.find("LFHCAL_7") != std::string::npos)
    return 147;

  else if (projname.find("BARR_0") != std::string::npos)
    return 150;
  else if (projname.find("BARR_1") != std::string::npos)
    return 151;
  else if (projname.find("BARR_2") != std::string::npos)
    return 152;
  else if (projname.find("BARR_3") != std::string::npos)
    return 153;
  else if (projname.find("BARR_4") != std::string::npos)
    return 154;
  else if (projname.find("SVTX_0") != std::string::npos)
    return 155;
  else if (projname.find("SVTX_1") != std::string::npos)
    return 156;
  else if (projname.find("SVTX_2") != std::string::npos)
    return 157;
  else if (projname.find("SVTX_3") != std::string::npos)
    return 158;
  else if (projname.find("SVTX_4") != std::string::npos)
    return 159;

  else if (projname.find("BST_0") != std::string::npos)
    return 160;
  else if (projname.find("BST_1") != std::string::npos)
    return 161;
  else if (projname.find("BST_2") != std::string::npos)
    return 162;
  else if (projname.find("BST_3") != std::string::npos)
    return 163;
  else if (projname.find("BST_4") != std::string::npos)
    return 164;

  else if (projname.find("BARR") != std::string::npos)
    return 170;
  else if (projname.find("SVTX") != std::string::npos)
    return 171;
  else
    return -1;
  return -1;
}

std::string EventEvaluatorEIC::GetProjectionNameFromIndex(int projindex)
{
  switch (projindex)
  {
  case 0:
    return "FTTL_0";
  case 1:
    return "FTTL_1";
  case 2:
    return "FTTL_2";
  case 3:
    return "ETTL_0";
  case 4:
    return "ETTL_1";
  case 5:
    return "LFHCAL";
  case 6:
    return "FEMC";
  case 7:
    return "CTTL_0";
  case 8:
    return "CTTL_1";

  case 35:
    return "hpDIRC";
  case 36:
    return "mRICH_Plane";
  case 37:
    return "dRICH_Plane";

  case 50:
    return "FST_0";
  case 51:
    return "FST_1";
  case 52:
    return "FST_2";
  case 53:
    return "FST_3";
  case 54:
    return "FST_4";
  case 55:
    return "FST_5";
  case 56 :
    return "FST_6";

  case 60:
    return "EHCAL";
  case 61:
    return "EEMC";
  case 62:
    return "HCALIN";
  case 63:
    return "HCALOUT";
  case 66:
    return "BECAL";

  case 70:
    return "ZDCsurrogate";
  case 71:
    return "rpTruth";
  case 72:
    return "b0Truth";
  case 73:
    return "rpTruth2";
  case 74:
    return "offMomTruth";

  case 79:
    return "TrackingService";

  case 90:
    return "BH_1";
  case 91:
    return "BH_FORWARD_PLUS";
  case 92:
    return "BH_FORWARD_NEG";

  case 110:
    return "RWELL_0";
  case 111:
    return "RWELL_1";
  case 112:
    return "RWELL_2";

  case 113:
    return "EST_0";
  case 114:
    return "EST_1";
  case 115:
    return "EST_2";
  case 116:
    return "EST_3";
  case 117:
    return "EST_4";
  case 118:
    return "EST_5";


  case 120:
    return "FGEM_0";
  case 121:
    return "FGEM_1";

  case 130:
    return "EGEM_0";
  case 131:
    return "EGEM_1";


  case 140:
    return "LFHCAL_0";
  case 141:
    return "LFHCAL_1";
  case 142:
    return "LFHCAL_2";
  case 143:
    return "LFHCAL_3";
  case 144:
    return "LFHCAL_4";
  case 145:
    return "LFHCAL_5";
  case 146:
    return "LFHCAL_6";
  case 147:
    return "LFHCAL_7";

  case 150:
    return "BARR_0";
  case 151:
    return "BARR_1";
  case 152:
    return "BARR_2";
  case 153:
    return "BARR_3";
  case 154:
    return "BARR_4";
  case 155:
    return "SVTX_0";
  case 156:
    return "SVTX_1";
  case 157:
    return "SVTX_2";
  case 158:
    return "SVTX_3";
  case 159:
    return "SVTX_4";

  case 160:
    return "BST_0";
  case 161:
    return "BST_1";
  case 162:
    return "BST_2";
  case 163:
    return "BST_3";
  case 164:
    return "BST_4";

  case 170:
    return "BARR";
  case 171:
    return "SVTX";
  default:
    return "NOTHING";
  }
}

int EventEvaluatorEIC::GetExponentFromProjectionIndex(int projindex)
{
  switch (projindex)
  {
  case 160: //"BST_0";
  case 155: //"SVTX_0";
    return 0;
  case 161: //"BST_1";
  case 156: //"SVTX_1";
    return 1;
  case 162: //"BST_2";
  case 157: //"SVTX_2";
    return 2;
  case 163: //"BST_3";
  case 150: //"BARR_0";
    return 3;
  case 164: // "BST_4";
  case 151: // "BARR_1";
    return 4;
  case 50: // "FST_0";
    return 5;
  case 51: // "FST_1";
    return 6;
  case 52: // "FST_2";
    return 7;
  case 53: // "FST_3";
    return 8;
  case 54: // "FST_4";
    return 9;
  case 55: // "FST_5";
    return 10;
  case 113: //"EST_0";
    return 11;
  case 114: //"EST_1";
    return 12;
  case 115: // "EST_2";
    return 13;
  case 116: // "EST_3";
    return 14;
  case 117: // "EST_4";
    return 15;
  case 110: //"RWELL_0";
    return 16;
  case 111: //"RWELL_1";
    return 17;
  case 112: //"RWELL_2";
    return 18;
  case 7: //"CTTL_0";
    return 19;
  case 0: //"FTTL_0";
    return 20;    
  case 1: // "FTTL_1";
    return 21;
  case 2: //  "FTTL_1";
    return 22;    
  case 3: // "ETTL_0";
    return 23;
  case 4: //"ETTL_1";
    return 24;

  default:
    return -1;
  }
}


void EventEvaluatorEIC::resetGeometryArrays()
{
  for (Int_t igeo = 0; igeo < _calo_towers_N; igeo++)
  {
    _calo_towers_iEta[_calo_towers_N] = -10000;
    _calo_towers_iPhi[_calo_towers_N] = -10000;
    _calo_towers_iL[_calo_towers_N] = -10000;
    _calo_towers_Eta[_calo_towers_N] = -10000;
    _calo_towers_Phi[_calo_towers_N] = -10000;
    _calo_towers_x[_calo_towers_N] = -10000;
    _calo_towers_y[_calo_towers_N] = -10000;
    _calo_towers_z[_calo_towers_N] = -10000;
  }
  _calo_ID = -1;
  _calo_towers_N = 0;
}
void EventEvaluatorEIC::resetBuffer()
{
  if (_do_store_event_info)
  {
    _cross_section = 0;
    _event_weight = 0;
    _n_generator_accepted = 0;
    if (Verbosity() > 0)
      std::cout << "\t... event info variables reset" << endl;
  }
  if (_do_VERTEX)
  {
    _vertex_x = -1000;
    _vertex_y = -1000;
    _vertex_z = -1000;
    _vertex_NCont = 0;
    _vertex_true_x = -1000;
    _vertex_true_y = -1000;
    _vertex_true_z = -1000;
    if (Verbosity() > 0)
      std::cout << "\t... vertex variables reset" << endl;
  }
  if (_do_HITS)
  {
    _nHitsLayers = 0;
    for (Int_t ihit = 0; ihit < _maxNHits; ihit++)
    {
      _hits_layerID[ihit] = 0;
      _hits_trueID[ihit] = 0;
      _hits_x[ihit] = 0;
      _hits_y[ihit] = 0;
      _hits_z[ihit] = 0;
      _hits_x2[ihit] = 0;
      _hits_y2[ihit] = 0;
      _hits_z2[ihit] = 0;
      _hits_t[ihit] = 0;
      _hits_edep[ihit] = 0;
      _hits_lightyield[ihit] = 0;
      _hits_isAbsorber[ihit] = 0;
    }
    if (Verbosity() > 0)
      std::cout << "\t... hit variables reset" << endl;
  }
  if (_do_BECAL)
  {
    _nTowers_BECAL = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_BECAL_E[itow] = 0;
      _tower_BECAL_iEta[itow] = 0;
      _tower_BECAL_iPhi[itow] = 0;
      _tower_BECAL_trueID[itow] = 0;
    }
    if (Verbosity() > 0)
      std::cout << "\t... BECAL variables reset" << endl;
  }
  if (_do_FEMC)
  {
    _nTowers_FEMC = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_FEMC_E[itow] = 0;
      _tower_FEMC_iEta[itow] = 0;
      _tower_FEMC_iPhi[itow] = 0;
      _tower_FEMC_trueID[itow] = 0;
    }
    if (_do_CLUSTERS)
    {
      _nclusters_FEMC = 0;
      for (Int_t itow = 0; itow < _maxNclusters; itow++)
      {
        _cluster_FEMC_E[itow] = 0;
        _cluster_FEMC_Eta[itow] = 0;
        _cluster_FEMC_Phi[itow] = 0;
        _cluster_FEMC_NTower[itow] = 0;
        _cluster_FEMC_trueID[itow] = 0;
      }
    }
    if (Verbosity() > 0)
      std::cout << "\t... FEMC variables reset" << endl;
  }
  if (_do_HCALIN)
  {
    _nTowers_HCALIN = 0;
    for (Int_t itow = 0; itow < _maxNTowersCentral; itow++)
    {
      _tower_HCALIN_E[itow] = 0;
      _tower_HCALIN_iEta[itow] = 0;
      _tower_HCALIN_iPhi[itow] = 0;
      _tower_HCALIN_trueID[itow] = 0;
    }
    if (_do_CLUSTERS)
    {
      _nclusters_HCALIN = 0;
      for (Int_t itow = 0; itow < _maxNclusters; itow++)
      {
        _cluster_HCALIN_E[itow] = 0;
        _cluster_HCALIN_Eta[itow] = 0;
        _cluster_HCALIN_Phi[itow] = 0;
        _cluster_HCALIN_NTower[itow] = 0;
        _cluster_HCALIN_trueID[itow] = 0;
      }
    }
    if (Verbosity() > 0)
      std::cout << "\t... HCALIN variables reset" << endl;
  }
  if (_do_HCALOUT)
  {
    if (Verbosity() > 0)
      std::cout << "\t... resetting HCALOUT variables" << endl;

    _nTowers_HCALOUT = 0;
    for (Int_t itow = 0; itow < _maxNTowersCentral; itow++)
    {
      _tower_HCALOUT_E[itow] = 0;
      _tower_HCALOUT_iEta[itow] = 0;
      _tower_HCALOUT_iPhi[itow] = 0;
      _tower_HCALOUT_trueID[itow] = 0;
    }
    if (_do_CLUSTERS)
    {
      _nclusters_HCALOUT = 0;
      for (Int_t itow = 0; itow < _maxNclusters; itow++)
      {
        _cluster_HCALOUT_E[itow] = 0;
        _cluster_HCALOUT_Eta[itow] = 0;
        _cluster_HCALOUT_Phi[itow] = 0;
        _cluster_HCALOUT_NTower[itow] = 0;
        _cluster_HCALOUT_trueID[itow] = 0;
      }
    }
    if (Verbosity() > 0)
      std::cout << "\t... HCALOUT variables reset" << endl;
  }

  if (_do_EEMC)
  {
    if (Verbosity() > 0)
      std::cout << "\t... resetting EEMC variables" << endl;
  
    _nTowers_EEMC = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_EEMC_E[itow] = 0;
      _tower_EEMC_iEta[itow] = 0;
      _tower_EEMC_iPhi[itow] = 0;
      _tower_EEMC_trueID[itow] = 0;
    }
    if (_do_CLUSTERS)
    {
      _nclusters_EEMC = 0;
      for (Int_t itow = 0; itow < _maxNclusters; itow++)
      {
        _cluster_EEMC_E[itow] = 0;
        _cluster_EEMC_Eta[itow] = 0;
        _cluster_EEMC_Phi[itow] = 0;
        _cluster_EEMC_NTower[itow] = 0;
        _cluster_EEMC_trueID[itow] = 0;
      }
    }
    if (Verbosity() > 0)
      std::cout << "\t... EEMC variables reset" << endl;  
  }

  if (_do_DRCALO)
  {
    if (Verbosity() > 0)
      std::cout << "\t... resetting DRCALO variables" << endl;

    _nTowers_DRCALO = 0;
    for (Int_t itow = 0; itow < _maxNTowersDR; itow++)
    {
      _tower_DRCALO_E[itow] = 0;
      _tower_DRCALO_NScint[itow] = 0;
      _tower_DRCALO_NCerenkov[itow] = 0;
      _tower_DRCALO_iEta[itow] = 0;
      _tower_DRCALO_iPhi[itow] = 0;
      _tower_DRCALO_trueID[itow] = 0;
    }
    if (Verbosity() > 0)
      std::cout << "\t... DRCALO variables reset" << endl;
  }

  if (_do_LFHCAL)
  {
    if (Verbosity() > 0)
      std::cout << "\t... resetting LFHCAL variables" << endl;
  
    _nTowers_LFHCAL = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_LFHCAL_E[itow] = 0;
      _tower_LFHCAL_iEta[itow] = 0;
      _tower_LFHCAL_iPhi[itow] = 0;
      _tower_LFHCAL_iL[itow] = 0;
      _tower_LFHCAL_trueID[itow] = 0;
    }
    if (Verbosity() > 0)
      std::cout << "\t... LFHCAL variables reset" << endl;
  }
  
  if (_do_TRACKS)
  {
    if (Verbosity() > 0)
      std::cout << "\t... resetting Track variables" << endl;
  
    _nTracks = 0;
    for (Int_t itrk = 0; itrk < _maxNTracks; itrk++)
    {
      _track_ID[itrk] = 0;
      _track_charge[itrk] = 0;
      _track_trueID[itrk] = 0;
      _track_px[itrk] = 0;
      _track_py[itrk] = 0;
      _track_pz[itrk] = 0;
      _track_x[itrk] = 0;
      _track_y[itrk] = 0;
      _track_z[itrk] = 0;
      _track_ndf[itrk] = 0;
      _track_chi2[itrk] = 0;
      _track_nHits[itrk]        = 0;
      _track_hitsEncoded[itrk]  = 0;
      _track_dca[itrk] = 0;
      _track_dca_2d[itrk] = 0;
      _track_source[itrk] = 0;

      // Default to not performing PID matching. Set default LL to same value for all PIDs
      _track_pion_LL[itrk] = -100;
      _track_kaon_LL[itrk] = -100;
      _track_proton_LL[itrk] = -100;
    }
    if (_do_PROJECTIONS)
    {
      _nProjections = 0;
      for (Int_t iproj = 0; iproj < _maxNProjections; iproj++)
      {
        _track_ProjLayer[iproj] = -1;
        _track_ProjTrackID[iproj] = 0;
        _track_TLP_x[iproj] = 0;
        _track_TLP_y[iproj] = 0;
        _track_TLP_z[iproj] = 0;
        _track_TLP_t[iproj] = 0;
        _track_TLP_px[iproj] = 0;
        _track_TLP_py[iproj] = 0;
        _track_TLP_pz[iproj] = 0;
        _track_TLP_true_x[iproj] = 0;
        _track_TLP_true_y[iproj] = 0;
        _track_TLP_true_z[iproj] = 0;
        _track_TLP_true_t[iproj] = 0;
      }
    }
    if (Verbosity() > 0)
      std::cout << "\t... track variables reset" << endl;
  }
  if (_do_MCPARTICLES)
  {
    _nMCPart = 0;
    for (Int_t imcpart = 0; imcpart < _maxNMCPart; imcpart++)
    {
      _mcpart_ID[imcpart] = 0;
      _mcpart_ID_parent[imcpart] = 0;
      _mcpart_PDG[imcpart] = 0;
      _mcpart_E[imcpart] = 0;
      _mcpart_px[imcpart] = 0;
      _mcpart_py[imcpart] = 0;
      _mcpart_pz[imcpart] = 0;
      _mcpart_x[imcpart] = 0;
      _mcpart_y[imcpart] = 0;
      _mcpart_z[imcpart] = 0;
      _mcpart_BCID[imcpart] = -10;
    }
  }

  if (_do_HEPMC)
  {
    _nHepmcp = 0;
    _hepmcp_vtx_x = 0;
    _hepmcp_vtx_y = 0;
    _hepmcp_vtx_z = 0;
    _hepmcp_vtx_t = 0;
    for (Int_t iHepmcp = 0; iHepmcp < _maxNHepmcp; iHepmcp++)
    {
      _hepmcp_E[iHepmcp] = 0;
      _hepmcp_PDG[iHepmcp] = 0;
      _hepmcp_px[iHepmcp] = 0;
      _hepmcp_py[iHepmcp] = 0;
      _hepmcp_pz[iHepmcp] = 0;
      _hepmcp_status[iHepmcp] = -10;
      _hepmcp_BCID[iHepmcp] = 0;
      _hepmcp_m2[iHepmcp] = 0;
      _hepmcp_m1[iHepmcp] = 0;
    }
    if (Verbosity() > 0)
      std::cout << "\t... MC variables reset" << endl;
  }
}
