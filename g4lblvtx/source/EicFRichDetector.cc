//____________________________________________________________________________..
//
// This is a working template for the G4 Construct() method which needs to be implemented
// We wedge a method between the G4 Construct() to enable volume hierarchies on the macro
// so here it is called ConstructMe() but there is no functional difference
// Currently this installs a simple G4Box solid, creates a logical volume from it
// and places it. Put your own detector in place (just make sure all active volumes
// get inserted into the m_PhysicalVolumesSet)
// 
// Rather than using hardcoded values you should consider using the parameter class
// Parameter names and defaults are set in EicFRichSubsystem::SetDefaultParameters()
// Only parameters defined there can be used (also to override in the macro)
// to avoids typos.
// IMPORTANT: parameters have no inherent units, there is a convention (cm/deg)
// but in any case you need to multiply them here with the correct CLHEP/G4 unit 
// 
// The place where you put your own detector is marked with
// //begin implement your own here://
// //end implement your own here://
// Do not forget to include the G4 includes for your volumes
//____________________________________________________________________________..

#include "EicFRichDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>

class G4VSolid;
class PHCompositeNode;

using namespace std;
double bp_r(double z_cm){return 0.05025461*z_cm-0.180808;}
double rmax(double z_cm){return 0.6624*z_cm;}
//G4Material * element_material( string identifier );
//void addDetectorSection( G4LogicalVolume *logicWorld , string name , double z_pos , double thick , string material , string color);
//____________________________________________________________________________..
EicFRichDetector::EicFRichDetector(PHG4Subsystem *subsys,
		PHCompositeNode *Node,
		PHParameters *parameters,
		const std::string &dnam)
	: PHG4Detector(subsys, Node, dnam)
	  , m_Params(parameters)
{
}

//_______________________________________________________________
int EicFRichDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
	set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if (iter != m_PhysicalVolumesSet.end())
	{
		return 1;
	}
	return 0;
}

//_______________________________________________________________
void EicFRichDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
	//begin implement your own here://
	// Do not forget to multiply the parameters with their respective CLHEP/G4 unit !

	// This class describes the EIC HM RICH detector from information provided by Evaristo Cisbani <evaristo.cisbani@roma1.infn.it>
	// The length in z for the gas volume is projected to be between 100 - 150 cm

	// ---------------- // ---------------- // ---------------- // ---------------- // ---------------- // ---------------- // ---------------- //
	// Parameters
	double overall_z_pos = 150 * cm;
	double overall_gas_length = 144 * cm;
	// ----------------
	// Mylar entrance
	double z_m_1 = overall_z_pos;
	double t_m_1 = 0.02*cm;
	// ----------------
	// Aerogel
	double z_aero = z_m_1+t_m_1;
	double t_aero = 4.*cm;
	// ----------------
	// PMMA
	double z_PMMA = z_aero+t_aero; 
	double t_PMMA = 2*mm;
	// ----------------
	// C2F6
	double z_C2F6 = z_PMMA+t_PMMA;
	double t_C2F6 = overall_gas_length;
	// ----------------
	// Mirror layer 1
	double z_mr_1 = z_C2F6+t_C2F6;
	double t_mr_1 = 0.1*mm;
	// ----------------
	// Mirror layer 2
	double z_mr_2 = z_mr_1+t_mr_1;
	double t_mr_2 = 0.05*mm;
	// ----------------
	// Mylar exit
	double z_m_2 = z_mr_2+t_mr_2;
	double t_m_2 = 0.02*cm;
	// ---------------- // ---------------- // ---------------- // ---------------- // ---------------- // ---------------- // ---------------- //
	addDetectorSection( logicWorld , "RICH_mylar_ent" , z_m_1  , t_m_1  , "mylar"   , "cyan"    );
	addDetectorSection( logicWorld , "RICH_aerogel"   , z_aero , t_aero , "aerogel" , "gray"    );
	addDetectorSection( logicWorld , "RICH_PMMA"      , z_PMMA , t_PMMA , "PMMA"    , "green"   );
	addDetectorSection( logicWorld , "RICH_C2F6"      , z_C2F6 , t_C2F6 , "C2F6"    , "magenta" );
	addDetectorSection( logicWorld , "RICH_mirror_l1" , z_mr_1 , t_mr_1 , "SiO2"    , "white"   );
	addDetectorSection( logicWorld , "RICH_mirror_l2" , z_mr_2 , t_mr_2 , "cf_epo"  , "yellow"  );
	addDetectorSection( logicWorld , "RICH_mylar_ext" , z_m_2  , t_m_2  , "mylar"   , "cyan"    );

	//end implement your own here://
	return;
}
// ======================================================================================================
void EicFRichDetector::Print(const std::string &what) const
{
	cout << "EicFRich Detector:" << endl;
	if (what == "ALL" || what == "VOLUME")
	{
		cout << "Version 0.1" << endl;
		cout << "Parameters:" << endl;
		m_Params->Print();
	}
	return;
}
// ======================================================================================================
G4Material *  EicFRichDetector::element_material( std::string identifier ){
	G4Material * G4_mat = G4Material::GetMaterial("G4_AIR");
	G4double density;
	G4int ncomponents, natoms;

	if(identifier=="mylar"){
		G4_mat = G4Material::GetMaterial("G4_MYLAR");
	}
	else if(identifier=="C2F6"){	
		G4_mat = new G4Material("C2F6", density = 0.0057 * g / cm3, ncomponents = 2);
		G4_mat->AddElement(G4Element::GetElement("C"), natoms = 2);
		G4_mat->AddElement(G4Element::GetElement("F"), natoms = 6);
	}
	else if(identifier=="PMMA"){
		G4_mat = new G4Material("PMMA", density = 1.18 * g / cm3, ncomponents = 3);
		G4_mat->AddElement(G4Element::GetElement("C"), 3.6 / (3.6 + 5.7 + 1.4));
		G4_mat->AddElement(G4Element::GetElement("H"), 5.7 / (3.6 + 5.7 + 1.4));
		G4_mat->AddElement(G4Element::GetElement("O"), 1.4 / (3.6 + 5.7 + 1.4));
	}
	else if(identifier=="SiO2"){
		G4_mat = new G4Material("SiO2", density = 2.5 * g / cm3 , ncomponents = 2);
		G4_mat->AddElement(G4Element::GetElement("Si"), natoms = 1);
		G4_mat->AddElement(G4Element::GetElement("O" ), natoms = 2);
	}
	else if(identifier=="aerogel"){
		G4Material *SiO2Aerogel = new G4Material("aerogel_SiO2", density = 2.5 * g / cm3 , ncomponents = 2);
		SiO2Aerogel->AddElement(G4Element::GetElement("Si"), 1);
		SiO2Aerogel->AddElement(G4Element::GetElement("O" ), 2);

		G4Material *air = G4Material::GetMaterial("G4_AIR");

		G4double fracMass;
		G4_mat = new G4Material("aerogel", density = 0.094 * g / cm3 , ncomponents = 2);
		G4_mat->AddMaterial(air        , fracMass = 96.*perCent);
		G4_mat->AddMaterial(SiO2Aerogel, fracMass =  4.*perCent);
	}
	else if(identifier=="cf_epo"){
		G4String symbol;
		G4Element* elH  = new G4Element("Hydrogen",symbol="H" , 1., 1.01*g/mole);
		G4Element* elC  = new G4Element("Carbon"  ,symbol="C" , 6., 12.01*g/mole);
		G4Element* elN  = new G4Element("Nitrogen",symbol="N" , 7., 14.01*g/mole);
		G4Element* elO  = new G4Element("Oxygen"  ,symbol="O" , 8., 16.00*g/mole);

		G4Material *Epoxy = new G4Material("Epoxy",  density = 1.16*g/cm3, natoms=4);
		Epoxy->AddElement(elH, 32); // Hydrogen
		Epoxy->AddElement(elN,  2); // Nitrogen
		Epoxy->AddElement(elO,  4); // Oxygen
		Epoxy->AddElement(elC, 15); // Carbon

		G4_mat = new G4Material("CarbonFiber",  density =  1.750*g/cm3, natoms=2);
		G4_mat->AddMaterial(G4Material::GetMaterial("G4_C"), 74.5*perCent);  // Carbon
		G4_mat->AddMaterial(Epoxy,                           25.5*perCent);  // Epoxy (scotch)
	}
	return G4_mat;
}
// ======================================================================================================
void  EicFRichDetector::addDetectorSection( G4LogicalVolume *logicWorld , std::string name , double z_pos , double thick , std::string material , std::string color){

	double z_det   [] = {z_pos,z_pos+thick};
	double rin [2] = {0};
	double rout[2] = {0};

	for(int i = 0 ; i < 2 ; i++){
		rin [i] = bp_r(z_det[i]);
		rout[i] = rmax(z_det[i]);
	}
	G4Material * G4_mat = element_material( material );

	G4RotationMatrix *rotm = new G4RotationMatrix();
	rotm->rotateX(0);
	rotm->rotateY(0);
	rotm->rotateZ(0);

	G4Color G4_color = G4Color(G4Colour::White());
	if     (color=="cyan"   ) G4_color = G4Color(G4Colour::Cyan());
	else if(color=="yellow" ) G4_color = G4Color(G4Colour::Yellow());
	else if(color=="green"  ) G4_color = G4Color(G4Colour::Green());
	else if(color=="gray"   ) G4_color = G4Color(G4Colour::Gray());
	else if(color=="magenta") G4_color = G4Color(G4Colour::Magenta());

	G4VSolid *G4_polycone = new G4Polycone(name,0,360*degree,2,z_det,rin,rout);
	G4LogicalVolume *logical = new G4LogicalVolume(G4_polycone,G4_mat, "EicFRichLogical");
	G4VisAttributes *vis_1 = new G4VisAttributes(G4_color);
	vis_1->SetForceSolid(true);
	logical->SetVisAttributes(vis_1);

	G4VPhysicalVolume *phy_1 = new G4PVPlacement(rotm, G4ThreeVector(0,0,0), logical , "EicFRich", logicWorld, 0, false, OverlapCheck());

	// add it to the list of placed volumes so the IsInDetector method picks them up
	m_PhysicalVolumesSet.insert(phy_1);
}
