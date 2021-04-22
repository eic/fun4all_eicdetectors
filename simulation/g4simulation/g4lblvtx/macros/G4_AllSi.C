#ifndef MACRO_G4ALLSI_C
#define MACRO_G4ALLSI_C

#include <g4lblvtx/AllSiliconTrackerSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4lblvtx/SimpleNtuple.h>
#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)
// ======================================================================================================================
void load_AllSi_geom(PHG4Reco *g4Reco,int det_ver=3){
	// Loading All-Si Tracker from dgml file
	AllSiliconTrackerSubsystem *allsili = new AllSiliconTrackerSubsystem();
	allsili->set_string_param("GDMPath",Form("detector/genfitGeom_AllSi_v%i.gdml",det_ver));

	if(det_ver==1){

		allsili->AddAssemblyVolume("BEAMPIPE"); // Load beampipe from the gdml file
	}
	else{
		PipeInit(); // Load beampipe from Fun4All rather than from gdml file
		double pipe_radius = 0;
		pipe_radius = Pipe(g4Reco,pipe_radius);
	}

	allsili->AddAssemblyVolume("VST");      // Barrel
	allsili->AddAssemblyVolume("FST");      // Forward disks
	allsili->AddAssemblyVolume("BST");      // Backward disks

	// this is for plotting single logical volumes for debugging and geantino scanning they end up at the center, you can plot multiple
	// but they end up on top of each other. They cannot coexist with the assembly volumes, the code will quit if you try to use both.
	// allsili->AddLogicalVolume("BstContainerVolume04");
	// allsili->AddLogicalVolume("FstContainerVolume00");
	// allsili->AddLogicalVolume("FstChipAssembly37");
	// allsili->AddLogicalVolume("VstStave00");

	allsili->SuperDetector("LBLVTX");
	allsili->SetActive();          // this saves hits in the MimosaCore volumes
	allsili->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(allsili);
}
// ======================================================================================================================
void add_AllSi_to_kalman(PHG4TrackFastSim *kalman, double pixel_size = 10., int det_ver=3){
	double um_to_cm = 1E-04; // Conversion factor from um to cm
	char nodename[100];
	int nBarrel = 6;
	int nDisks = 5;

	if(det_ver==3){
		nBarrel = 7;
		nDisks = 7;
	}
	
	for (int i=10; i<10+nBarrel; i++){ // CENTRAL BARREL
		sprintf(nodename,"G4HIT_LBLVTX_CENTRAL_%d", i);
		kalman->add_phg4hits(
				nodename,                               // const std::string& phg4hitsNames
				PHG4TrackFastSim::Cylinder,             // const DETECTOR_TYPE phg4dettype
				999.,                                   // radial-resolution [cm] (this number is not used in cylindrical geometry)
				pixel_size*um_to_cm/sqrt(12.),          // azimuthal (arc-length) resolution [cm]
				pixel_size*um_to_cm/sqrt(12.),          // longitudinal (z) resolution [cm]
				1,                                      // efficiency (fraction)
				0                                       // hit noise
				);
	}
	for (int i=20; i<20+nDisks; i++){ // FORWARD DISKS
		sprintf(nodename,"G4HIT_LBLVTX_FORWARD_%d", i);
		kalman->add_phg4hits(
				nodename,                               // const std::string& phg4hitsNames
				PHG4TrackFastSim::Vertical_Plane,       // const DETECTOR_TYPE phg4dettype
				pixel_size*um_to_cm/sqrt(12.),          // radial-resolution [cm]
				pixel_size*um_to_cm/sqrt(12.),          // azimuthal (arc-length) resolution [cm]
				999.,                                   // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
				1,                                      // efficiency (fraction)
				0                                       // hit noise
				);
	}
	for (int i=30; i<30+nDisks; i++){ // BACKWARD DISKS
		sprintf(nodename,"G4HIT_LBLVTX_BACKWARD_%d", i);
		kalman->add_phg4hits(
				nodename,                               // const std::string& phg4hitsNames
				PHG4TrackFastSim::Vertical_Plane,       // const DETECTOR_TYPE phg4dettype
				pixel_size*um_to_cm/sqrt(12.),          // radial-resolution [cm]
				pixel_size*um_to_cm/sqrt(12.),          // azimuthal (arc-length) resolution [cm]
				999.,                                   // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
				1,                                      // efficiency (fraction)
				0                                       // hit noise
				);
	}
}
// ======================================================================================================================
void add_AllSi_hits(SimpleNtuple * hits, int det_ver=3){
	char nodename[100];
	int nBarrel = 6;
        int nDisks = 5;
       
        if(det_ver==3){
                nBarrel = 7;
                nDisks = 7;
        }

        //  hits->AddNode("ABSORBER_LBLVTX",0); // hits in the passive volumes
        for (int i = 10; i < 10+nBarrel; i++){sprintf(nodename, "LBLVTX_CENTRAL_%d", i);        hits->AddNode(nodename, i);} // hits in the  MimosaCore volumes
        for (int i = 20; i < 20+nDisks ; i++){sprintf(nodename, "LBLVTX_FORWARD_%d", i);        hits->AddNode(nodename, i);} // hits in the  MimosaCore volumes
        for (int i = 30; i < 30+nDisks ; i++){sprintf(nodename, "LBLVTX_BACKWARD_%d",i);        hits->AddNode(nodename, i);} // hits in the  MimosaCore volumes
}
#endif
