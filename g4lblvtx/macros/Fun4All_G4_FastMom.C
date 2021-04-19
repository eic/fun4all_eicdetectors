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
#include <phpythia8/PHPythia8.h>
#include <g4histos/G4HitNtuple.h>
#include <g4lblvtx/AllSiliconTrackerSubsystem.h>
#include <g4main/HepMCNodeReader.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <phool/recoConsts.h>

#include "G4_Jets.C"
#include "G4_Bbc.C"
#include "G4_Global.C"
#include "G4_Pipe_EIC.C"
#include "G4_AllSi.C"

#include <g4lblvtx/G4LBLVtxSubsystem.h>
#include <g4lblvtx/SimpleNtuple.h>
#include <g4lblvtx/TrackFastSimEval.h>
#include <myjetanalysis/MyJetAnalysis_AllSi.h>
#include <g4lblvtx/PHG4ParticleGenerator_flat_pT.h>
R__LOAD_LIBRARY(libmyjetanalysis.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)
R__LOAD_LIBRARY(libPHPythia8.so)

void Fun4All_G4_FastMom(
		int nEvents = -1,			// number of events
		const char *outputFile = "out_allSi",	// output filename
		const char *genpar = "pi-",		// particle species to simulate with the simple generators
		const int det_ver = 2,			// version of detector geometry
		const double pixel_size = 10.)		// pixel length (um)
{
	// ======================================================================================================
	// Input from the user
	const int particle_gen = 5;	// 1 = particle generator, 2 = particle gun, 3 = simple event generator, 4 = pythia8 e+p collision, 5 = particle generator flat in pT
	const int magnetic_field = 4;	// 1 = uniform 1.5T, 2 = uniform 3.0T, 3 = sPHENIX 1.4T map, 4 = Beast 3.0T map
	bool DISPLACED_VERTEX = false;	// this option exclude vertex in the track fitting and use RAVE to reconstruct primary and 2ndary vertexes
	bool do_projections = false;	// Project momentum vectors to surfaces defined below
	// ======================================================================================================
	// Parameters for projections
	string projname1   = "DIRC";            // Cylindrical surface object name
	double projradius1 = 50.;               // [cm] 
	double length1     = 400.;              // [cm]
	// ---
	double thinness    = 0.1;               // black hole thickness, needs to be taken into account for the z positions
	// ---
	string projname2   = "FOR";             // Forward plane object name
	double projzpos2   = 130+thinness/2.;   // [cm]
	double projradius2 = 50.;               // [cm]
	// ---
	string projname3   = "BACK";            // Backward plane object name
	double projzpos3   = -(130+thinness/2.);// [cm]
	double projradius3 = 50.;               // [cm]
	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	// If you want to fix the random seed for reproducibility
	// recoConsts *rc = recoConsts::instance();
	// rc->set_IntFlag("RANDOMSEED", 12345);
	// ======================================================================================================
	// Particle Generation
	if(particle_gen<4) cout << "Particle that will be generated: " << std::string(genpar) << endl;
	// --------------------------------------------------------------------------------------
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string(genpar));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);                    // Vertex generation range
	gen->set_mom_range(.00001,30.);		// Momentum generation range in GeV/c
	gen->set_z_range(0.,0.);
	gen->set_eta_range(-4,4);		// Detector coverage corresponds to |η|< 4
	gen->set_phi_range(0.,2.*TMath::Pi());
	// --------------------------------------------------------------------------------------
	// Particle generator flat in pT
	PHG4ParticleGenerator_flat_pT *gen_pT = new PHG4ParticleGenerator_flat_pT();
	gen_pT->set_name(std::string(genpar));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen_pT->set_vtx(0,0,0);                    // Vertex generation range
	gen_pT->set_pT_range(.00001,30.);         // Momentum generation range in GeV/c
	gen_pT->set_z_range(0.,0.);
	gen_pT->set_eta_range(-4,4);               // Detector coverage corresponds to |η|< 4
	gen_pT->set_phi_range(0.,2.*TMath::Pi());
	// --------------------------------------------------------------------------------------
	// Particle Gun Setup
	PHG4ParticleGun *gun = new PHG4ParticleGun();
	gun->set_name(std::string(genpar));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ...
	gun->set_vtx(0,0,0);
	gun->set_mom(0,1,0);
	// --------------------------------------------------------------------------------------
	// Simple event generator
	PHG4SimpleEventGenerator *segen = new PHG4SimpleEventGenerator();
	segen->add_particles(std::string(genpar),100); // 100 pion option
	segen->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,PHG4SimpleEventGenerator::Uniform,PHG4SimpleEventGenerator::Uniform);
	segen->set_vertex_distribution_mean(0.0, 0.0, 0.0);
	segen->set_vertex_distribution_width(0.0, 0.0, 5.0);
	segen->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
	segen->set_vertex_size_parameters(0.0, 0.0);
	segen->set_eta_range(-4.,4.);
	segen->set_phi_range(-TMath::Pi(),TMath::Pi());
	segen->set_p_range(.00001,30.);
	//segen->Embed(2);
	segen->Verbosity(0);
	// --------------------------------------------------------------------------------------
	// Pythia 8 events
	bool do_pythia8_jets = false;
	if     (particle_gen==1){se->registerSubsystem(  gen); cout << "Using particle generator"     << endl;}
	else if(particle_gen==2){se->registerSubsystem(  gun); cout << "Using particle gun"           << endl;}
	else if(particle_gen==3){se->registerSubsystem(segen); cout << "Using simple event generator" << endl;}
	else if(particle_gen==4){
		do_pythia8_jets = true;
		gSystem->Load("libPHPythia8.so");

		PHPythia8 *pythia8 = new PHPythia8(); // see coresoftware/generators/PHPythia8 for example config
		pythia8->set_config_file("phpythia8.cfg"); // example configure files : https://github.com/sPHENIX-Collaboration/coresoftware/tree/master/generators/PHPythia8 
		se->registerSubsystem(pythia8);

		// read-in HepMC events to Geant4 if there are any pythia8 produces hepmc records, so this is needed to read the above generated pythia8 events
		HepMCNodeReader *hr = new HepMCNodeReader();
		se->registerSubsystem(hr);
	}
	else if (particle_gen==5){se->registerSubsystem(gen_pT); cout << "Using particle generator flat in pT"  << endl;}
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
	// Physics list (default list is "QGSP_BERT")
	//g4Reco->SetPhysicsList("FTFP_BERT_HP"); // This list is slower and only useful for hadronic showers.
	// ======================================================================================================
	load_AllSi_geom(g4Reco, det_ver);	// Loading All-Si Tracker and beampipe geometries

	// ======================================================================================================
	if(do_projections){
		PHG4CylinderSubsystem *cyl;
		cyl = new PHG4CylinderSubsystem(projname1,0);
		cyl->set_double_param("length", length1);
		cyl->set_double_param("radius", projradius1); // dirc radius
		cyl->set_double_param("thickness", 0.1); // needs some thickness
		cyl->set_string_param("material", "G4_AIR");
		cyl->SetActive(1);
		cyl->SuperDetector(projname1);
		cyl->BlackHole();
		cyl->set_color(1,0,0,0.7); //reddish
		g4Reco->registerSubsystem(cyl);

		cyl = new PHG4CylinderSubsystem(projname2,0);
		cyl->set_double_param("length", thinness);
		cyl->set_double_param("radius", 2); // beampipe needs to fit here
		cyl->set_double_param("thickness", projradius2); // 
		cyl->set_string_param("material", "G4_AIR");
		cyl->set_double_param("place_z", projzpos2);
		cyl->SetActive(1);
		cyl->SuperDetector(projname2);
		cyl->BlackHole();
		cyl->set_color(0,1,1,0.3); //reddish
		g4Reco->registerSubsystem(cyl);

		cyl = new PHG4CylinderSubsystem(projname3,0);
		cyl->set_double_param("length", thinness);
		cyl->set_double_param("radius", 2); // beampipe needs to fit here
		cyl->set_double_param("thickness", projradius3); // 
		cyl->set_string_param("material", "G4_AIR");
		cyl->set_double_param("place_z", projzpos3);
		cyl->SetActive(1);
		cyl->SuperDetector(projname3);
		cyl->BlackHole();
		cyl->set_color(0,1,1,0.3); //reddish
		g4Reco->registerSubsystem(cyl);
	}
	// ======================================================================================================

	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem(g4Reco);

	// ======================================================================================================
	if(do_pythia8_jets){
		Bbc_Reco();
		Global_Reco();
	}

	// ======================================================================================================
	// Fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	double um_to_cm = 1E-04; // Conversion factor from um to cm
	char nodename[100];
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("SVTX");
	kalman->set_trackmap_out_name("SvtxTrackMap");

	add_AllSi_to_kalman( kalman , pixel_size , det_ver );	// Add All-Silicon tracker to Kalman filter

	// Projections  
	if(do_projections){
		kalman->add_cylinder_state(projname1, projradius1);     // projection on cylinder (DIRC)
		kalman->add_zplane_state  (projname2, projzpos2  );     // projection on vertical planes
		kalman->add_zplane_state  (projname3, projzpos3  );     // projection on vertical planes
	}

	if(DISPLACED_VERTEX){
		// use very loose vertex constraint (1cm in sigma) to allow reco of displaced vertex
		kalman->set_use_vertex_in_fitting(true);
		kalman->set_vertex_xy_resolution(1);
		kalman->set_vertex_z_resolution(1);
		kalman->enable_vertexing(true);
		kalman->set_vertex_min_ndf(2); // Only tracks with number of degrees of freedom greater than this value are used to fit the vertex
	}
	/*
	   else{
	   kalman->set_use_vertex_in_fitting(false);
	   kalman->set_vertex_xy_resolution(0);
	   kalman->set_vertex_z_resolution(0);
	   kalman->enable_vertexing(false); // this is false by default
	   kalman->set_vertex_min_ndf(2);
	   }
	*/
	//kalman -> Verbosity(10);
	se->registerSubsystem(kalman);
	// -----------------------------------------------------
	// INFO: The resolution numbers above correspond to:
	// 20e-4/sqrt(12) cm = 5.8e-4 cm, to simulate 20x20 um

	// ======================================================================================================
	TrackFastSimEval *fast_sim_eval = new TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(TString(outputFile)+B_label+"_FastTrackingEval.root");
	if(do_projections){
		fast_sim_eval->AddProjection(projname1);
		fast_sim_eval->AddProjection(projname2);
		fast_sim_eval->AddProjection(projname3);
	}
	se->registerSubsystem(fast_sim_eval);

	// ======================================================================================================
	// resonstruct jets after the tracking
	if(do_pythia8_jets) Jet_Reco();

	SimpleNtuple *hits = new SimpleNtuple("Hits");
	add_AllSi_hits(hits,det_ver);	// Add All-Silicon tracker hits
	se->registerSubsystem(hits);

	// ======================================================================================================
	if(do_pythia8_jets)
		Jet_Eval(string(outputFile) + "_g4jet_eval.root");
	// ======================================================================================================

	///////////////////////////////////////////
	// IOManagers...
	///////////////////////////////////////////
	const std::string dst_name = std::string(outputFile)+"_G4LBLVtx.root";
	//Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",TString(outputFile)+"_G4LBLVtx.root");
	Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",dst_name);
	out->Verbosity(0);
	se->registerOutputManager(out);

	Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
	se->registerInputManager(in);

	if(do_pythia8_jets){
		gSystem->Load("libmyjetanalysis");
		std::string jetoutputFile = std::string(outputFile) + std::string("_electrons+jets.root");
		MyJetAnalysis_AllSi *myJetAnalysis = new MyJetAnalysis_AllSi("AntiKt_Track_r10","AntiKt_Truth_r10",jetoutputFile.data());	
		//MyJetAnalysis_AllSi *myJetAnalysis = new MyJetAnalysis_AllSi("AntiKt_Track_r04","AntiKt_Truth_r04",jetoutputFile.data());
		//MyJetAnalysis_AllSi *myJetAnalysis = new MyJetAnalysis_AllSi("AntiKt_Track_r02","AntiKt_Truth_r02",jetoutputFile.data());
		se->registerSubsystem(myJetAnalysis);
	}

	if (nEvents <= 0) return;

	se->run(nEvents);
	se->End();
	delete se;

	gSystem->Exit(0);
}
