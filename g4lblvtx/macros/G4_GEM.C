#ifndef MACRO_G4GEM_C
#define MACRO_G4GEM_C

#include "GlobalVariables.C"

#include <g4detectors/PHG4SectorSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

	int make_GEM_station(string name, PHG4Reco *g4Reco, double zpos, double etamin, double etamax, const int N_Sector = 8);
	void AddLayers_MiniTPCDrift(PHG4SectorSubsystem *gem);

	namespace Enable
{
	bool EGEM = true;
	bool FGEM = true;
}  // namespace Enable
// ======================================================================================================================
void EGEM_Init(){
	BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 80.);
	// extends only to -z
	BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -160.);
}
// ======================================================================================================================
void FGEM_Init(){
	BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 150.);
	BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 282.);
}
// ======================================================================================================================
void EGEMSetup(PHG4Reco *g4Reco){
	/* Careful with dimensions! If GEM station volumes overlap, e.g. with TPC volume, they will be drawn in event display but will NOT register any hits.
	 * Geometric constraints:
	 * TPC length = 211 cm --> from z = -105.5 to z = +105.5
	 */
	float thickness = 3.;	
	make_GEM_station("EGEM", g4Reco, -180.0 , -1.15 , -3.9);
}
// ======================================================================================================================
void FGEMSetup(PHG4Reco *g4Reco){
	float thickness = 3.;
	make_GEM_station("FGEM", g4Reco, 300.0 , 1.4, 3.69);
}
// ======================================================================================================================
void addPassiveMaterial(PHG4Reco *g4Reco){
	float z_pos = 130.0;

	// This is a mockup calorimeter in the forward (hadron-going) direction
	PHG4CylinderSubsystem *cyl_f = new PHG4CylinderSubsystem("CALO_FORWARD_PASSIVE",0);
	cyl_f->set_double_param("length", 5);		// Length in z direction in cm
	cyl_f->set_double_param("radius", z_pos*0.0503-0.180808); // beampipe needs to fit here
	cyl_f->set_double_param("thickness", 43); // 
	cyl_f->set_string_param("material", "G4_Al");
	cyl_f->set_double_param("place_z", z_pos);
	//cyl_f->SetActive(1);
	cyl_f->SuperDetector("passive_F");
	//cyl_f->set_color(0,1,1,0.3); //reddish
	g4Reco->registerSubsystem(cyl_f);

	
	// This is a mockup calorimeter in the backward (electron-going) direction
	PHG4CylinderSubsystem * cyl_b = new PHG4CylinderSubsystem("CALO_BACKWARD_PASSIVE",0);
	cyl_b->set_double_param("length", 5);	// Length in z direction in cm
	cyl_b->set_double_param("radius",abs(-z_pos*0.030-0.806));	// beampipe needs to fit here
	cyl_b->set_double_param("thickness", 43); // 
	cyl_b->set_string_param("material", "G4_Al");
	cyl_b->set_double_param("place_z", -z_pos);
	//cyl_b->SetActive(1);
	cyl_b->SuperDetector("passive_B");
	//cyl_b->set_color(0,1,1,0.3); //reddish
	g4Reco->registerSubsystem(cyl_b);
	
}
// ======================================================================================================================
int make_GEM_station(string name, PHG4Reco *g4Reco, double zpos, double etamin,
		double etamax, const int N_Sector = 8)
{
	double polar_angle = 0;

	if (zpos < 0){
		zpos = -zpos;
		polar_angle = M_PI;
	}
	if (etamax < etamin){
		double t = etamax;
		etamax = etamin;
		etamin = t;
	}

	PHG4SectorSubsystem *gem;
	gem = new PHG4SectorSubsystem(name);

	gem->SuperDetector(name);

	gem->get_geometry().set_normal_polar_angle(polar_angle);
	gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
	gem->get_geometry().set_min_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamax));
	gem->get_geometry().set_max_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamin));
	gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
	gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
	gem->get_geometry().set_N_Sector(N_Sector);
	gem->get_geometry().set_material("G4_METHANE");
	gem->OverlapCheck(Enable::OVERLAPCHECK);

	AddLayers_MiniTPCDrift(gem);
	gem->get_geometry().AddLayers_HBD_GEM();
	g4Reco->registerSubsystem(gem);
	return 0;
}
// ======================================================================================================================
//! Add drift layers to mini TPC
void AddLayers_MiniTPCDrift(PHG4SectorSubsystem *gem)
{
	assert(gem);

	const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
	const double mm = 0.1 * cm;
	const double um = 1e-3 * mm;

	//  const int N_Layers = 70; // used for mini-drift TPC timing digitalization
	const int N_Layers = 1;  // simplified setup
	const double thickness = 2 * cm;

	gem->get_geometry().AddLayer("EntranceWindow", "G4_MYLAR", 25 * um, false, 100);
	gem->get_geometry().AddLayer("Cathode", "G4_GRAPHITE", 10 * um, false, 100);

	for (int d = 1; d <= N_Layers; d++)
	{
		ostringstream s;
		s << "DriftLayer_";
		s << d;

		gem->get_geometry().AddLayer(s.str(), "G4_METHANE", thickness / N_Layers, true);
	}
}
#endif
