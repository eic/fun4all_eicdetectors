/*
================================================================================================================
The purpose of this code is to have a version of the all-silicon tracker that is simplified so that we can study
different variations of the geometry quickly. Specifically, I wrote this code to study the impact that changing
the material budget of different regions of the detector would have on different resolutions.
================================================================================================================
*/
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
#include <g4lblvtx/PHG4ParticleGenerator_flat_pT.h>
#include <g4lblvtx/AllSi_Al_support_Subsystem.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

void Fun4All_G4_simplified_v2(
			int nEvents = -1,			// number of events
			double vtx_matBud  = 0.05,		// % X/X0 (material budget of vertexing layers)
			double barr_matBud = 0.55,		// % X/X0 (material budget of middle layers)
			double disk_matBud = 0.25,		// % X/X0 (material budget of disk layers)
			double pmin = 0., 			// GeV/c
			double pmax = 30., 			// GeV/c
			int magnetic_field = 4, 		// Magnetic field setting
			TString out_name = "out_vtx_study")	// output filename
{	
	// ======================================================================================================
	// Input from the user
	const int particle_gen = 5;     // 1 = particle generator, 2 = particle gun, 3 = simple event generator, 4 = pythia8 e+p collision, 5 = particle generator flat in pT
	double pix_size_vtx = 10.; // um - size of pixels in vertexing layers
	double pix_size_bar = 10.; // um - size of pixels in barrel layers
	double pix_size_dis = 10.; // um - size of pixels in disk layers
	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	// If you want to fix the random seed for reproducibility
	// recoConsts *rc = recoConsts::instance();
	// rc->set_IntFlag("RANDOMSEED", 12345);
	// ======================================================================================================
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string("pi-"));	// geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);			// Vertex generation range
	gen->set_mom_range(pmin,pmax);		// Momentum generation range in GeV/c
	gen->set_z_range(0.,0.);
	gen->set_eta_range(0.,4.);
	gen->set_phi_range(0,2.*TMath::Pi());
	// --------------------------------------------------------------------------------------
	// Particle generator flat in pT
	PHG4ParticleGenerator_flat_pT *gen_pT = new PHG4ParticleGenerator_flat_pT();
	gen_pT->set_name(std::string("pi-"));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen_pT->set_vtx(0,0,0);                    // Vertex generation range
	gen_pT->set_pT_range(.00001,5.);         // Momentum generation range in GeV/c
	gen_pT->set_z_range(0.,0.);
	gen_pT->set_eta_range(-4,4);               // Detector coverage corresponds to |η|< 4
	gen_pT->set_phi_range(0.,2.*TMath::Pi());
	// ======================================================================================================
	if     (particle_gen==1){se->registerSubsystem(  gen); cout << "Using particle generator"     << endl;}
	else if(particle_gen==5){se->registerSubsystem(gen_pT); cout << "Using particle generator flat in pT"  << endl;}
	else{ cout << "Particle generator option requested has not been implemented. Bailing out!" << endl; exit(0); }
	// ======================================================================================================
	PHG4Reco *g4Reco = new PHG4Reco();
	//g4Reco->SetWorldMaterial("G4_Galactic");	
	// ======================================================================================================
	// Magnetic field setting
	TString B_label;
	if(magnetic_field==1){          // uniform 1.5T
		B_label = "_B_1.5T";
		g4Reco->set_field(1.5);
	}
	else if(magnetic_field==2){     // uniform 3.0T
		B_label = "_B_3.0T";
		g4Reco->set_field(3.0);
	}
	else if(magnetic_field==3){     // sPHENIX 1.4T map
		B_label = "_sPHENIX";
		g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"), PHFieldConfig::kField2D);
		g4Reco->set_field_rescale(-1.4/1.5);
	}
	else if(magnetic_field==4){     // Beast 3.0T map
		B_label = "_Beast";
		g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/mfield.4col.dat"), PHFieldConfig::kFieldBeast);
	}
	else{                           // The user did not provide a valid B field setting
		cout << "User did not provide a valid magnetic field setting. Set 'magnetic_field'. Bailing out!" << endl;
	}	
	// ======================================================================================================
	// Detector setup
	PHG4CylinderSubsystem *cyl;
	//---------------------------
	// Vertexing
	double si_vtx_r_pos[] = {3.30,5.70};
	const int nVtxLayers = sizeof(si_vtx_r_pos)/sizeof(*si_vtx_r_pos);
	double si_z_vtxlength[] = {30.,30.};
	double si_thick_vtx = vtx_matBud/100.*9.37;

	for (int ilayer = 0; ilayer < nVtxLayers ; ilayer++){
		cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
		cyl->set_string_param("material" , "G4_Si"               );
		cyl->set_double_param("radius"   , si_vtx_r_pos[ilayer]  );
		cyl->set_double_param("thickness", si_thick_vtx          );
		cyl->set_double_param("place_z"  , 0                     );
		cyl->set_double_param("length"   , si_z_vtxlength[ilayer]);
		cyl->SetActive();
		cyl->SuperDetector("SVTX");
		g4Reco->registerSubsystem(cyl);
	}
	//---------------------------
	// Barrel
	double si_r_pos[] = {21.,22.68,39.3,43.23};
	const int nTrckLayers = sizeof(si_r_pos)/sizeof(*si_r_pos);
	double si_z_length[] = {54.,60.,105.,114.};
	double si_thick_bar = barr_matBud/100.*9.37;

	for (int ilayer = 0; ilayer < nTrckLayers ; ilayer++){
		cyl = new PHG4CylinderSubsystem("BARR", ilayer);
		cyl->set_string_param("material" , "G4_Si"            );
		cyl->set_double_param("radius"   , si_r_pos[ilayer]   );
		cyl->set_double_param("thickness", si_thick_bar       );
		cyl->set_double_param("place_z"  , 0                  );
		cyl->set_double_param("length"   , si_z_length[ilayer]);
		cyl->SetActive();
		cyl->SuperDetector("BARR");
		cyl->set_color(0,0.5,1);
		g4Reco->registerSubsystem(cyl);	
	}
	//---------------------------
	// Disks
	double si_z_pos[] = {-121.,-97.,-73.,-49.,-25.,25.,49.,73.,97.,121.};
	double si_r_max[10] = {0};
	double si_r_min[10] = {0};
	double si_thick_disk = disk_matBud/100.*9.37;
	for(int i = 0 ; i < 10 ; i++){
		si_r_max[i] = TMath::Min(43.23,18.5*abs(si_z_pos[i])/si_z_pos[5]);

		if(si_z_pos[i]>66.8&&si_z_pos[i]>0) si_r_min[i] = (0.05025461*si_z_pos[i]-0.180808);
		else if(si_z_pos[i]>0) si_r_min[i] = 3.18;
		else if(si_z_pos[i]<-79.8&&si_z_pos[i]<0) si_r_min[i] = (-0.0297039*si_z_pos[i]+0.8058281);
		else si_r_min[i] = 3.18;

		si_r_max[i] -= si_r_min[i];
	}

	for (int ilayer = 0; ilayer < 10; ilayer++){
		cyl = new PHG4CylinderSubsystem("FBVS", ilayer);
		cyl->set_string_param("material" , "G4_Si"         );
		cyl->set_double_param("radius"   , si_r_min[ilayer]);
		cyl->set_double_param("thickness", si_r_max[ilayer]);
		cyl->set_double_param("place_z"  , si_z_pos[ilayer]);
		cyl->set_double_param("length"   , si_thick_disk   );
		cyl->SetActive();
		cyl->SuperDetector("FBST");
		cyl->set_color(1,0,0);
		g4Reco->registerSubsystem(cyl);
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

	// ------------
	// Al Support Structure
	AllSi_Al_support_Subsystem *Al_supp = new AllSi_Al_support_Subsystem("Al_supp");
	g4Reco->registerSubsystem(Al_supp);	
	// ------------	

	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem(g4Reco);

	//---------------------------
	// fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	//---------------------------
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("BARR");
	kalman->set_trackmap_out_name("SvtxTrackMap");

	// add Vertexing Layers
	kalman->add_phg4hits(
			"G4HIT_SVTX",				// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,					// radial-resolution [cm]
			pix_size_vtx/10000./sqrt(12.),		// azimuthal-resolution [cm]
			pix_size_vtx/10000./sqrt(12.),		// z-resolution [cm]
			1,					// efficiency,
			0					// noise hits
			);

	// add Barrel Layers
	kalman->add_phg4hits(
			"G4HIT_BARR",                   	// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,                           	// radial-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      	// azimuthal-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      	// z-resolution [cm]
			1,                              	// efficiency,
			0                               	// noise hits
			);

	//  add Disk Layers
	kalman->add_phg4hits(
			"G4HIT_FBST",				// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Vertical_Plane,
			pix_size_dis/10000./sqrt(12.),		// radial-resolution [cm]
			pix_size_dis/10000./sqrt(12.),		// azimuthal-resolution [cm]
			999.,                       		// z-resolution [cm]
			1,                          		// efficiency,
			0                           		// noise hits
			);	

	//kalman->Verbosity(10);
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_vertex_xy_resolution(0);
	kalman->set_vertex_z_resolution(0);
	kalman->enable_vertexing(false); // this is false by default
	kalman->set_vertex_min_ndf(2);

	se->registerSubsystem(kalman);

	std::string outputFile = (std::string)(out_name)+std::string(B_label)+"_FastSimEval.root";

	PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(outputFile);
	se->registerSubsystem(fast_sim_eval);

	// ======================================================================================================
	// IOManagers...
	const std::string dst_name = std::string(out_name)+std::string(B_label)+"_G4LBLVtx.root";	
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
