#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <phool/PHRandomSeed.h>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/HepMCNodeReader.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <phool/recoConsts.h>
#include "G4Setup_sPHENIX.C"
#include "G4_Global.C"
#include "G4_DSTReader.C"
#include "DisplayOn.C"

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libphhepmc.so)

#endif

using namespace std;

int Fun4All_G4_sPHENIX_particle_gen(
			const int nEvents = 1,
			const char *outputFile = "G4sPHENIX.root",
			const char *genpar = "pi-")
{
	const bool usegun = true;
	bool do_pipe = true;
	bool do_tracking = true;
	bool do_global_fastsim = false;

	//---------------
	// Load libraries
	//---------------
	gSystem->Load("libfun4all.so");
	gSystem->Load("libg4detectors.so");	
	gSystem->Load("libg4testbench.so");
	gSystem->Load("libg4eval.so");
	gSystem->Load("libg4intt.so");

	// establish the geometry and reconstruction setup
	gROOT->LoadMacro("G4Setup_sPHENIX.C");
	G4Init(do_tracking,false,false,false,false,false,do_pipe,false,false);

	int absorberactive = 1;  // set to 1 to make all absorbers active volumes

	//  const string magfield = "1.5"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
	const string magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"); // default map from the calibration database
	const float magfield_rescale = -1.4 / 1.5;                                     // scale the map to a 1.4 T field

	//---------------
	// Fun4All server
	//---------------
	bool display_on = false;
	if(display_on) gROOT->LoadMacro("DisplayOn.C");

	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(0);

	// just if we set some flags somewhere in this macro
	recoConsts *rc = recoConsts::instance();
	//  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
	//  rc->set_IntFlag("RANDOMSEED", 12345);

	//-----------------
	// Event generation
	//-----------------
	if (usegun)
	{
		PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
		gen->set_name(std::string(genpar));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
		gen->set_vtx(0,0,0);                    // Vertex generation range
		gen->set_mom_range(.00001,8.);         // Momentum generation range in GeV/c
		gen->set_z_range(0.,0.);
		gen->set_eta_range(-1.5,1.5);               // Detector coverage corresponds to |η|< 4
		gen->set_phi_range(0.,2.*TMath::Pi());
		se->registerSubsystem(gen);
	}

	//---------------------
	// Detector description
	//---------------------
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
	G4Setup(absorberactive, magfield, EDecayType::kAll,
			do_tracking, false , false , false , false , false , do_pipe, false , false , magfield_rescale);
#else
	G4Setup(absorberactive, magfield, TPythia6Decayer::kAll,
			do_tracking, false , false , false , false , false , do_pipe, false , false , magfield_rescale);
#endif
	Tracking_Cells();
	Tracking_Clus();
	Tracking_Reco();

	//-----------------
	// Global Vertexing
	//-----------------
	if (do_global_fastsim){
		gROOT->LoadMacro("G4_Global.C");
		Global_FastSim();
	}

	Tracking_Eval(string(outputFile) + "_g4svtx_eval.root");

	//--------------
	// IO management
	//--------------
	Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
	se->registerInputManager(in);

	//-----------------
	// Event processing
	//-----------------
	if (nEvents <= 0) return 0;

	if(display_on)
	{
		DisplayOn();
		// prevent macro from finishing so can see display
		int i;
		cout << "***** Enter any integer to proceed" << endl;
		cin >> i;
	}

	se->run(nEvents);

	//-----
	// Exit
	//-----
	se->End();
	std::cout << "All done" << std::endl;
	delete se;
	gSystem->Exit(0);
	return 0;
}
