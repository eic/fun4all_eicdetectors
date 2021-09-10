#include "PHG4TRDDetector.h"
#include <phparameter/PHParameters.h>
#include <g4main/PHG4Detector.h>                // for PHG4Detector

#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4PVDivision.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4Element.hh>
#include <Geant4/G4VPhysicalVolume.hh>          // for G4VPhysicalVolume
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>              // for G4ThreeVector
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>                             // for operator<<, endl, bas...

class PHCompositeNode;
class PHG4Subsystem;

using namespace std;


PHG4TRDDetector::PHG4TRDDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , Phys(nullptr)
  , fPhysicsRadiator(nullptr)   
  , TRD_det_Phys(nullptr)
  , m_Active(m_Params->get_int_param("active"))
  , m_AbsorberActive(m_Params->get_int_param("absorberactive"))
  , m_Layer(lyr)
{

}

//_______________________________________________________________
//_______________________________________________________________
int PHG4TRDDetector::IsInTRD(const G4VPhysicalVolume *volume) const
{
  if(m_Active)
    {
      if (volume == fPhysicsRadiator )
	{
	  return -1;
	}
    }
  
  if(m_AbsorberActive)
    {
      if(volume == Gas_Active)
	{
	  return 1;
	}
    }
  return 0;

}

//_______________________________________________________________
void PHG4TRDDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  // Generic mother volume 
  //All other components of TRD will be in it 
  G4Material *TRDMaterial = G4Material::GetMaterial(m_Params->get_string_param("material"));
  
  double RIn = m_Params->get_double_param("RIn")*cm; // 20.* cm
  double ROut = m_Params->get_double_param("ROut")*cm; // 200 * cm
  double ThicknessZ = m_Params->get_double_param("ThicknessZ")*cm; // 30. *cm
  double PosZ = m_Params->get_double_param("PosZ")*cm; // none defined
  
  cout << "Mother CENTER Pos  ( z ) :" << " " << PosZ << " Mother length " << ThicknessZ << endl;
  
  G4Tubs* Solid = new G4Tubs("TRD_GVol_Solid", RIn, ROut, ThicknessZ / 2., 0., 360 * deg);
  G4LogicalVolume *Logic = new G4LogicalVolume(Solid, TRDMaterial, "TRD_GVol_Logic",0 ,0,0);
  G4VisAttributes *attr_TRD_GVol = new G4VisAttributes(G4Color(0.3, 0.5, 0.9, 0.9));
  //attr_TRD_GVol->SetColor(G4Color::Green());
  attr_TRD_GVol->SetForceSolid(true);
  Logic->SetVisAttributes(attr_TRD_GVol);
  Phys =  new G4PVPlacement(0, G4ThreeVector(0, 0, PosZ), Logic, "H_CAP_TRD_Physics", logicWorld, 0, false, OverlapCheck());
  

  //Radiator construction ---------
  double det_RIn = m_Params->get_double_param("det_RIn")*cm;
  double det_ROut = m_Params->get_double_param("det_ROut")*cm;
  //Phase space distribution of radiators -------------
  double fGasGap = 0.0600*cm ; // for ZEUS 300 publication
  double fDetGap = 0.001 * cm;
  double fRadThickness = 0.0020 * cm;    // 16 um // ZEUS NIMA 323 (1992) 135-139, D=20um, dens.= 0.1 g/cm3
  //double fRadThick = 18.*cm  - fGasGap + fDetGap ;
  double fRadThick = 10.*cm  - fGasGap + fDetGap ;
  // int fFoilNumber = 0;
  //fFoilNumber = fRadThick / (fRadThickness + fGasGap);
  
  double fRadZ ; // none defined
  //fRadZ = -ThicknessZ/2. + fRadThick/2. + 2.* cm;
  fRadZ = (-ThicknessZ/2. + fRadThick/2. );
 
  double  foilGasRatio = fRadThickness / (fRadThickness + fGasGap);

  // ------ Define materials for radiator----------
 
  //G4Material *Mylar = G4Material::GetMaterial("Mylar");
  G4Material *Air = G4Material::GetMaterial("G4_AIR");
  //G4Material *Al = G4Material::GetMaterial("G4_Al");
  //G4Material *CH2 = G4Material::GetMaterial("G4_CH2");
  //G4Material *He = G4Material::GetMaterial("He");
  
  double a_c = 12.0107*g/mole;
  G4Element *Carbon = new G4Element("Carbon", "C" ,6,  a_c) ;
  double a_h = 1.01*g/mole;
  G4Element *Hydrogen = new G4Element("Hydrogen", "H" ,1, a_h) ;
  
  double density_ch2 = 0.935 * g/cm3 ; // g/cm^3
  G4Material *CH2 = new G4Material("CH2", density_ch2, 2);
  CH2->AddElement(Carbon, 1);
  CH2->AddElement(Hydrogen, 2);
  

  double foilDensity = 0.91 * g / cm3;  // CH2 1.39*g/cm3; // Mylar //  0.534*g/cm3; //Li
  double gasDensity = 1.2928 * mg / cm3; // Air // 1.977*mg/cm3; // CO2 0.178*mg/cm3; // He
  double totDensity = foilDensity * foilGasRatio + gasDensity * (1.0 - foilGasRatio);
  double fractionFoil = foilDensity * foilGasRatio / totDensity;
  double fractionGas = gasDensity * (1.0 - foilGasRatio) / totDensity;

  
  G4Material *radiatorMat0 = new G4Material("radiatorMat0", totDensity, 2);
  radiatorMat0->AddMaterial(CH2, fractionFoil);
  radiatorMat0->AddMaterial(Air, fractionGas);
  double NewDensity = 0.083 * (g / cm3);
  G4Material *radiatorMat = new G4Material("radiatorMat", NewDensity, 1);
  radiatorMat->AddMaterial(radiatorMat0, 1.);
 
  //------------------------end --material -------------------------------


  //Define all sorts of G4 volumes 
  G4Tubs *fSolidRadiator = new G4Tubs("TRD_Radiator_Solid", det_RIn, det_ROut, 0.5 *fRadThick, 0., 360 * deg);
  G4LogicalVolume *fLogicRadiator = new G4LogicalVolume(fSolidRadiator, radiatorMat, "TRD_Radiator_Logic",0,0,0); 

  G4VisAttributes *TRD_rad_att = new G4VisAttributes();
  TRD_rad_att->SetColour(G4Colour::Blue());
  TRD_rad_att->SetVisibility(true);
  TRD_rad_att->SetForceSolid(true);
  fLogicRadiator->SetVisAttributes(TRD_rad_att);

   fPhysicsRadiator =  new G4PVPlacement(0, G4ThreeVector(0, 0, fRadZ),  fLogicRadiator,  "TRD_Radiator_Phys", Logic, 0 , false, OverlapCheck()); // Placing the Radiator in global detector volume "Phys". Phys has been placed in world volume earlier
   /*
  G4Region *fRadRegion = new G4Region("XTRradiator");
  fRadRegion->AddRootLogicalVolume(fLogicRadiator);
   */
  

   // Absorber (Xe gas with MPGD)
   
   //Cover for gas 
   //Mylar cover 
   double window_th =  0.002*cm;
   //double off = 0.001*cm; 
 // double window_pos_Z =  -1.0 * det_ThicknessZ/2 +  off;
   double window_pos_Z = fRadZ + fRadThick/2. + window_th/2. ;
   G4Material *window_Material = G4Material::GetMaterial("G4_MYLAR");
   
   G4Tubs *MPGD_win_Solid = new G4Tubs("MPGD_win_Solid", det_RIn, det_ROut, window_th/2. , 0. , 360 * deg);;
   G4LogicalVolume *MPGD_win_Logic =  new G4LogicalVolume(MPGD_win_Solid, window_Material, "MPGD_win_Logic",0, 0, 0);
   G4VisAttributes *MPGD_win_att = new G4VisAttributes();
   MPGD_win_att->SetColour(G4Colour::Gray());
   MPGD_win_att->SetVisibility(true);
   MPGD_win_att->SetForceSolid(true);
   MPGD_win_Logic->SetVisAttributes(MPGD_win_att);
   
   MPGD_win_Phys = new G4PVPlacement(0, G4ThreeVector(0,0, window_pos_Z), MPGD_win_Logic, "MPGD_win_Phys", Logic, 0, false, OverlapCheck()); 
  
        
     double det_ThicknessZ =2.5*cm;
   // double det_Pos_Z = fRadZ+ fRadThick/2. + det_ThicknessZ/2. ;
   double det_Pos_Z =  window_pos_Z + window_th/2. + det_ThicknessZ/2. ;
   //G4Material *det_Material = G4Material::GetMaterial("G4_Xe");
   G4Material *det_Material = G4Material::GetMaterial("G4_AIR");
   
   G4Tubs *TRD_det_Solid = new G4Tubs("TRD_det_Solid", det_RIn, det_ROut, det_ThicknessZ/2., 0., 360 * deg);;
   G4LogicalVolume *TRD_det_Logic =  new G4LogicalVolume(TRD_det_Solid, det_Material, "TRD_det_Logic",0, 0, 0);
   G4VisAttributes *TRD_det_att = new G4VisAttributes();
   TRD_det_att->SetColour(G4Colour::Red());
   TRD_det_att->SetVisibility(true);
   TRD_det_att->SetForceSolid(true);
   TRD_det_Logic->SetVisAttributes(TRD_det_att);
   
   TRD_det_Phys =  new G4PVPlacement(0, G4ThreeVector(0, 0, det_Pos_Z),  TRD_det_Logic,  "TRD_det_Phys", Logic, 0 , false, OverlapCheck()); // Placing the MPGD with Xe in global detector volume "Phys". Phys has been placed in world volume earlier
  
   cout  << " det RIN :" << det_RIn << " det ROUT :" << det_ROut <<  endl;
  
   
   //Define MPGD detector with drift cathode and R/o. These are the daughter of TRD_det_Logic
  

   //construct drift cathode

   double cat_th = 0.005*cm;
   double dead_ar = 0.01*cm;
   G4Material *G4_Al = G4Material::GetMaterial("G4_Al");
   G4Material *myCatMesh = new G4Material("myCatMesh",0.9*g/cm3, G4_Al, kStateSolid); // Density of Al reduced considering it's a mesh

   G4Tubs *Cathode = new G4Tubs("Cathode", det_RIn, det_ROut, cat_th/2., 0., 360 *deg);
   G4LogicalVolume *Cathode_Logic =  new G4LogicalVolume(Cathode, myCatMesh, " Cathode_Logic", 0, 0, 0);
   G4VisAttributes *Cathode_att = new G4VisAttributes();
   Cathode_att->SetColour(G4Colour::Magenta());
   Cathode_att->SetVisibility(true);
   Cathode_att->SetForceSolid(true);
   Cathode_Logic->SetVisAttributes(Cathode_att);
   double cat_pos  = -1.0*det_ThicknessZ/2. + dead_ar + cat_th/2.;
   
   Cathode_Phys =  new G4PVPlacement(0, G4ThreeVector(0, 0, cat_pos), Cathode_Logic, "Cathode_Phys", TRD_det_Logic, 0, false, OverlapCheck());

   //construct active gas volume upto top GEM
   double gas_thick =  2.0*cm;
   G4Material *gas_act = G4Material::GetMaterial("G4_Xe");
   G4Tubs *drift_gas = new G4Tubs("drift_gas", det_RIn, det_ROut, gas_thick/2., 0., 360*deg);
   G4LogicalVolume *gas_Logic = new G4LogicalVolume(drift_gas, gas_act, "gas_Logic", 0, 0, 0);
   G4VisAttributes *gas_att = new G4VisAttributes();
   gas_att->SetColour(G4Colour::Green());
   gas_att->SetVisibility(true);
   gas_att->SetForceSolid(true);
   gas_Logic->SetVisAttributes(gas_att);
   double gas_pos = cat_pos +cat_th/2.+gas_thick/2.; 

   Gas_Active = new G4PVPlacement(0, G4ThreeVector(0, 0, gas_pos), gas_Logic, "Gas_Active", TRD_det_Logic, 0, false, OverlapCheck());
   

  // construct GEMs and MMG and R/o
  //GEMs + top + below region dimensions
  double cu_th = 0.0005*cm;
  // double dr_gap = 2.0*cm;
  double tr_gap = 0.2*cm;
  double kap_th = 0.005*cm;
  G4Material *gem_cond = G4Material::GetMaterial("G4_Cu");
  G4Material *gem_diel = G4Material::GetMaterial("G4_KAPTON");
  
  //MMG + below dimensions
  double av_gap = 0.014*cm;
  double mesh_th =  0.0012*cm;

  G4Material *G4_Cr = G4Material::GetMaterial("G4_Cr");
  G4Material *G4_Fe = G4Material::GetMaterial("G4_Fe");
  G4Material *G4_Mn = G4Material::GetMaterial("G4_Mn");
  G4Material *G4_Ni = G4Material::GetMaterial("G4_Ni");
  G4Material *G4_Si = G4Material::GetMaterial("G4_Si");
  G4Material *G4_O = G4Material::GetMaterial("G4_O");
  G4Material *G4_H = G4Material::GetMaterial("G4_H");

  G4Material *myMMMesh = new G4Material( "myMMMesh", 2.8548*g/cm3, 5, kStateSolid);
  myMMMesh->AddMaterial( G4_Cr, 0.1900 );
  myMMMesh->AddMaterial( G4_Fe, 0.6800 );
  myMMMesh->AddMaterial( G4_Mn, 0.0200 );
  myMMMesh->AddMaterial( G4_Ni, 0.1000 );
  myMMMesh->AddMaterial( G4_Si, 0.0100 );

  //R/o dimensions
  double res_th = 0.0020* cm ;
  double pcb_th = 0.01*cm;
  double cu_st_th = 0.0012*cm;

  G4Material *G4_C = G4Material::GetMaterial("G4_C");
  G4Material *Reslay =  new G4Material("Reslay", 0.77906*g/cm3, G4_C, kStateSolid);
  
  G4Material *G4_Cu = G4Material::GetMaterial("G4_Cu");
  G4Material *MMstrips = new G4Material("MMstrips",  5.28414*g/cm3, G4_Cu, kStateSolid);

  
  G4Material *myFR4 = new G4Material("myFR4", 1.860*g/cm3, 4, kStateSolid);
  myFR4->AddMaterial( G4_C, 0.43550 );
  myFR4->AddMaterial( G4_H, 0.03650 );
  myFR4->AddMaterial( G4_O, 0.28120 );
  myFR4->AddMaterial( G4_Si, 0.24680 );
  
  G4Tubs *GEM_top_Solid = new G4Tubs("GEM_top_Solid", det_RIn, det_ROut, cu_th/2., 0., 360 * deg);
  G4LogicalVolume *GEM_top_Logic = new G4LogicalVolume(GEM_top_Solid,  gem_cond , "GEM_top_Logic", 0,  0, 0 );
  G4VisAttributes *GEM_top_att = new G4VisAttributes();
  GEM_top_att->SetColour(G4Colour::Yellow());
  GEM_top_att->SetVisibility(true);
  GEM_top_att->SetForceSolid(true);
  GEM_top_Logic->SetVisAttributes(GEM_top_att);
  //double gem_top_z = /*window_pos_Z + window_th/2. + dr_gap + cu_th/2.*/ -1.0*det_ThicknessZ/2. + dr_gap + cu_th/2. ; 
  //double gem_top_z = cat_pos + cat_th/2. + dr_gap + cu_th/2; 
  double gem_top_z = gas_pos + gas_thick/2. + cu_th/2; 
  GEM_top_Phys = new G4PVPlacement(0, G4ThreeVector(0, 0, gem_top_z), GEM_top_Logic, "GEM_top_Phys", TRD_det_Logic, 0, false, OverlapCheck());

  G4Tubs *GEM_diel_Solid = new G4Tubs("GEM_diel_Solid", det_RIn, det_ROut, kap_th/2., 0., 360 * deg);
  G4LogicalVolume *GEM_diel_Logic = new G4LogicalVolume(GEM_diel_Solid, gem_diel, "GEM_diel_Logic", 0, 0, 0 );
  G4VisAttributes *GEM_diel_att = new G4VisAttributes();
  GEM_diel_att->SetColour(G4Colour::Green());
  GEM_diel_att->SetVisibility(true);
  GEM_diel_att->SetForceSolid(true);
  GEM_diel_Logic->SetVisAttributes(GEM_diel_att);

  double diel_z = gem_top_z + cu_th/2. + kap_th/2. ; 

  GEM_diel_Phys = new G4PVPlacement(0, G4ThreeVector(0, 0, diel_z), GEM_diel_Logic, "GEM_diel_Phys", TRD_det_Logic, 0, false, OverlapCheck());

  G4Tubs *GEM_bottom_Solid = new G4Tubs("GEM_bottom_Solid", det_RIn, det_ROut, cu_th/2., 0., 360 * deg);
  G4LogicalVolume *GEM_bottom_Logic = new G4LogicalVolume(GEM_bottom_Solid, gem_cond, "GEM_bottom_Logic", 0, 0, 0 );
  G4VisAttributes *GEM_bottom_att = new G4VisAttributes();
  GEM_bottom_att->SetColour(G4Colour::Yellow());
  GEM_bottom_att->SetVisibility(true);
  GEM_bottom_att->SetForceSolid(true);
  GEM_bottom_Logic->SetVisAttributes(GEM_bottom_att);
  
  double gem_bot_z = diel_z + kap_th/2. + cu_th/2. ;
  GEM_bottom_Phys = new G4PVPlacement(0, G4ThreeVector(0, 0, gem_bot_z), GEM_bottom_Logic, "GEM_bottom_Phys", TRD_det_Logic, 0, false, OverlapCheck());

 
  G4Tubs *MMG_mesh_Solid = new G4Tubs("MMG_mesh_Solid", det_RIn, det_ROut, cu_th/2., 0., 360 * deg);
  G4LogicalVolume *MMG_mesh_Logic = new G4LogicalVolume(MMG_mesh_Solid, myMMMesh, "MMG_mesh_Logic", 0, 0, 0 );
  G4VisAttributes *MMG_mesh_att = new G4VisAttributes();
  MMG_mesh_att->SetColour(G4Colour::Brown());
  MMG_mesh_att->SetVisibility(true);
  MMG_mesh_att->SetForceSolid(true);
  MMG_mesh_Logic->SetVisAttributes(MMG_mesh_att);
  
  double mmg_mesh_z = gem_bot_z + cu_th/2. + tr_gap + mesh_th/2.;
  MMG_mesh_Phys = new G4PVPlacement(0, G4ThreeVector(0, 0, mmg_mesh_z), MMG_mesh_Logic, "MMG_mesh_Phys", TRD_det_Logic, 0, false, OverlapCheck());

  G4Tubs *Res_lay_Solid = new G4Tubs("Res_lay_Solid", det_RIn, det_ROut, res_th/2., 0., 360 * deg);
  G4LogicalVolume *Res_lay_Logic = new G4LogicalVolume(Res_lay_Solid,  Reslay , "Res_lay_Logic", 0,  0,  0 );
  G4VisAttributes *Res_lay_att = new G4VisAttributes();
  Res_lay_att->SetColour(G4Colour::Brown());
  Res_lay_att->SetVisibility(true);
  Res_lay_att->SetForceSolid(true);
  Res_lay_Logic->SetVisAttributes(Res_lay_att);
  
  double res_lay_z = mmg_mesh_z +  mesh_th/2. + av_gap + res_th/2. ;
  Res_lay_Phys = new G4PVPlacement(0, G4ThreeVector(0, 0, res_lay_z), Res_lay_Logic, "Res_lay_Phys", TRD_det_Logic, 0, false, OverlapCheck());

  G4Tubs *MMG_strips_Solid = new G4Tubs("MMG_strips_Solid", det_RIn, det_ROut,  cu_st_th/2., 0., 360. * deg);
  G4LogicalVolume *MMG_strips_Logic = new G4LogicalVolume(MMG_strips_Solid,  MMstrips, "MMG_strips_Logic", 0, 0, 0);
  G4VisAttributes *MMG_strips_att = new G4VisAttributes();
  MMG_strips_att->SetColour(G4Colour::Yellow());
  MMG_strips_att->SetVisibility(true);
  MMG_strips_att->SetForceSolid(true);
  MMG_strips_Logic->SetVisAttributes(MMG_strips_att);

  double mmg_str_z = res_lay_z + res_th/2. + cu_st_th/2. ; 
  MMG_strips_Phys = new G4PVPlacement(0, G4ThreeVector(0 , 0, mmg_str_z), MMG_strips_Logic, "MMG_strips_Phys", TRD_det_Logic, 0, false, OverlapCheck());
  
  G4Tubs *PCB_Solid = new G4Tubs("PCB_Solid", det_RIn, det_ROut,  pcb_th/2., 0., 360. * deg);
  G4LogicalVolume *PCB_Logic = new G4LogicalVolume(PCB_Solid,  myFR4, "PCB_Logic", 0, 0, 0);
  G4VisAttributes *PCB_att = new G4VisAttributes();
  PCB_att->SetColour(G4Colour::Green());
  PCB_att->SetVisibility(true);
  PCB_att->SetForceSolid(true);
  PCB_Logic->SetVisAttributes(PCB_att);

  double pcb_z = mmg_str_z + cu_st_th/2. +  pcb_th/2.; 
  PCB_Phys = new G4PVPlacement(0, G4ThreeVector(0 , 0, pcb_z), PCB_Logic, "PCB_Phys", TRD_det_Logic, 0, false, OverlapCheck());
  
  

  /*
    G4Region * fRegGasDet = new G4Region("XTRdEdxDetector"); 
    fRegGasDet->AddRootLogicalVolume(TRD_det_Logic);
  */
}
