#pragma once
#include <phgenfit/Track.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4histos/G4HitNtuple.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <phool/recoConsts.h>
#include <g4lblvtx/EicFRichSubsystem.h>

#include "G4_GEM.C"

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

void Fun4All_G4_simple_hadron_GEM(
			int nEvents = -1,		// number of events
			float phi_deg_min = 0,		// phi range
			TString out_name = "out_AKiselev_GEM")	// output filename
{
	float phi_min = phi_deg_min*TMath::Pi()/180.;
	TString phi_range = Form("_%.0f_%.0f",phi_deg_min,phi_deg_min+45.);
	TString outputFile = out_name+phi_range+"_FastSimEval.root";
	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	// ======================================================================================================
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string("pi-"));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);                    // Vertex generation range
	gen->set_mom_range(5.,50.);         // Momentum generation range in GeV/c
	gen->set_z_range(0.,0.);
	gen->set_eta_range(2,6);
	gen->set_phi_range(phi_min,phi_min+TMath::Pi()/4.);
	se->registerSubsystem(gen);
	// ======================================================================================================
	PHG4Reco *g4Reco = new PHG4Reco();
	g4Reco->SetWorldMaterial("G4_Galactic");	
	g4Reco->set_field(3.0);
	// ======================================================================================================
	double si_z_pos[] = {25.,49.,73.,97.,121};
	double si_max_r = 44.;
	double si_min_r = 0.0001;
	double si_thick = 0.02811;

	PHG4CylinderSubsystem *cyl;
	// here is our silicon:
	for (int ilayer = 0; ilayer < 5; ilayer++)
	{
		cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
		cyl->set_string_param("material" , "G4_Si"         );
		cyl->set_double_param("radius"   , si_min_r        );
		cyl->set_double_param("thickness", si_max_r        );
		cyl->set_double_param("place_z"  , si_z_pos[ilayer]);
		cyl->set_double_param("length"   , si_thick        );
		cyl->SetActive();
		cyl->SuperDetector("SVTX");
		g4Reco->registerSubsystem(cyl);
	}

	FGEM_Init();	// Loading forward GEM geometry
	FGEMSetup(g4Reco);
	// ------------
	// Forward RICH
	EicFRichSubsystem *RICH = new EicFRichSubsystem("RICH");
	g4Reco->registerSubsystem(RICH);

	//---------------------------
	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem(g4Reco);

	//---------------------------
	// fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	//---------------------------
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("SVTX");
	kalman->set_trackmap_out_name("SvtxTrackMap");

	//  add Si Trtacker
	kalman->add_phg4hits(
			"G4HIT_SVTX",               // const std::string& phg4hitsNames,
			PHG4TrackFastSim::Vertical_Plane,
			20./10000./sqrt(12.),       // radial-resolution [cm]
			20./10000./sqrt(12.),       // azimuthal-resolution [cm]
			999.,                       // z-resolution [cm]
			1,                          // efficiency,
			0                           // noise hits
			);
	
	kalman->add_phg4hits("G4HIT_FGEM",                  // const std::string& phg4hitsNames,
                        PHG4TrackFastSim::Vertical_Plane,   // const DETECTOR_TYPE phg4dettype,
                        50e-4,//1. / sqrt(12.),             // const float radres,
                        50e-4,                              // const float phires,
                        999.,                               // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
                        1,                                  // const float eff,
                        0                                   // const float noise
                        );
	
	//kalman->Verbosity(10);
	se->registerSubsystem(kalman);

	PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(outputFile);
	se->registerSubsystem(fast_sim_eval);

	// ======================================================================================================
	// IOManagers...
	const std::string dst_name = std::string(out_name+phi_range)+"_G4LBLVtx.root";	
	Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",dst_name);
	out->Verbosity(0);
	se->registerOutputManager(out);

	Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
	se->registerInputManager(in);

	if (nEvents <= 0) return;

	se->run(nEvents);
	se->End();
	delete se;

	gSystem->Exit(0);
}
