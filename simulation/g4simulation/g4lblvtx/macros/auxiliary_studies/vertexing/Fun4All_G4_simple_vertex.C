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

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

	void Fun4All_G4_simple_vertex(
			int nEvents = -1,		// number of events
			bool vtx_lyr_1 = true,
			bool vtx_lyr_2 = true,
			bool vtx_lyr_3 = true,
			double vtx_matBud = 0.05, //% X/X0
			double pix_size = 10.,
			TString out_name = "out_vtx_study")	// output filename
{	
	std::string outputFile = std::string(out_name)+"_FastSimEval.root";
	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	// ======================================================================================================
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string("pi-"));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);                    // Vertex generation range
	gen->set_mom_range(0.,8.);         // Momentum generation range in GeV/c
	gen->set_z_range(0.,0.);
	gen->set_eta_range(0,1);
	gen->set_phi_range(0,2*TMath::Pi());
	se->registerSubsystem(gen);
	// ======================================================================================================
	PHG4Reco *g4Reco = new PHG4Reco();
	g4Reco->SetWorldMaterial("G4_Galactic");	
	g4Reco->set_field(3.0);
	// ======================================================================================================
	// Detector setup
	double si_r_pos[] = {3.64,4.45,5.26,21.,22.68,39.3,43.23};
	const int nTrckLayers = sizeof(si_r_pos)/sizeof(*si_r_pos);
	double si_z_length[] = {14.,14.,14.,18.,20.,35.,38.};
	for(int i = 0 ; i < nTrckLayers ; i++) si_z_length[i] *= 3.;
	double si_thick_vtx = vtx_matBud/100.*9.37;
	double si_thick_bar = 0.55/100.*9.37;

	PHG4CylinderSubsystem *cyl;
	for (int ilayer = 0; ilayer < nTrckLayers ; ilayer++)
	{
		if(
				(ilayer==0&&vtx_lyr_1)||
				(ilayer==1&&vtx_lyr_2)||
				(ilayer==2&&vtx_lyr_3)||
				(ilayer>2)
		  ){
			cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
			cyl->set_string_param("material" , "G4_Si"         );
			cyl->set_double_param("radius"   , si_r_pos[ilayer]);

			if(ilayer<2)
				cyl->set_double_param("thickness", si_thick_vtx);
			else
				cyl->set_double_param("thickness", si_thick_bar);

			cyl->set_double_param("place_z"  , 0 );
			cyl->set_double_param("length"   , si_z_length[ilayer]    );
			cyl->SetActive();
			cyl->SuperDetector("SVTX");
			g4Reco->registerSubsystem(cyl);
		}
	}

	//---------------------------
	// mid-rapidity beryllium pipe
	double be_pipe_radius = 3.1000;
	double be_pipe_thickness = 3.1762 - be_pipe_radius;  // 760 um for sPHENIX
	double be_pipe_length_plus = 66.8;                   // +z beam pipe extend.
	double be_pipe_length_neg = -79.8;                   // -z beam pipe extend.
	double be_pipe_length = be_pipe_length_plus - be_pipe_length_neg;
	double be_pipe_center = 0.5 * (be_pipe_length_plus + be_pipe_length_neg);

	cyl = new PHG4CylinderSubsystem("BE_PIPE", 1);
	cyl->set_double_param("radius", be_pipe_radius);
	cyl->set_int_param("lengthviarapidity", 0);
	cyl->set_double_param("length", be_pipe_length);
	cyl->set_double_param("place_z", be_pipe_center);
	cyl->set_string_param("material", "G4_Be");
	cyl->set_double_param("thickness", be_pipe_thickness);
	cyl->SuperDetector("PIPE");
	g4Reco->registerSubsystem(cyl);
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
			"G4HIT_SVTX",			// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,				// radial-resolution [cm]
			pix_size/10000./sqrt(12.),	// azimuthal-resolution [cm]
			pix_size/10000./sqrt(12.),	// z-resolution [cm]
			1,				// efficiency,
			0				// noise hits
			);
	//kalman->Verbosity(10);
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_vertex_xy_resolution(0);
	kalman->set_vertex_z_resolution(0);
	kalman->enable_vertexing(false); // this is false by default
	kalman->set_vertex_min_ndf(2);

	se->registerSubsystem(kalman);

	PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(outputFile);
	se->registerSubsystem(fast_sim_eval);

	// ======================================================================================================
	// IOManagers...
	const std::string dst_name = std::string(out_name)+"_G4LBLVtx.root";	
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
