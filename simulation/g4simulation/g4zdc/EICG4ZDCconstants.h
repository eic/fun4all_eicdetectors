
#ifndef EICG4ZDCconstants_h
#define EICG4ZDCconstants_h 1

#include <Geant4/G4SystemOfUnits.hh>

//Crystal Towers
constexpr G4int nCTowerX=20;
constexpr G4int nCTowerY=20;
constexpr G4int nCTowerZ=1;
constexpr G4double CTower_X = 3. *cm;
constexpr G4double CTower_Y = 3. *cm;
constexpr G4double CTower_Z = 7. *cm;
constexpr G4double CTower_GAP = 3. *cm;

//MiniFoCal block
constexpr G4int NumberOfLayers = 22;
constexpr G4int NumberPIX = 2;
constexpr G4int NumberPAD = 20;
constexpr G4int NPadOnlyLayers=10;

//=================================================================
//PIX detector
constexpr G4int NpixX = 200;
constexpr G4int NpixY = 200;
constexpr G4double PIX_X = 3.*mm;   //changed from 0.1
constexpr G4double PIX_Y = 3.*mm;    //changed from 0.1
constexpr G4double PIX_Z = 0.30*mm;  //changed from 0.1
//Thicknesses
constexpr G4double PIX_Absorber_Z = 3.5*mm;
constexpr G4double PIX_Glue1_Z  = 0.11*mm; //between W-Si
constexpr G4double PIX_Glue2_Z  = 0.11*mm; //between Si+FPC
constexpr G4double PIX_FPC_Z    = 0.28*mm; //Readout
constexpr G4double PIX_AirGap   = 1.2*mm;
constexpr G4double PIX_Layer_Thickness = PIX_Absorber_Z + PIX_Glue1_Z + PIX_Z + PIX_Glue2_Z + PIX_FPC_Z + PIX_AirGap;

//=================================================================
//PAD detector
constexpr G4int NpadX = 60;
constexpr G4int NpadY = 60;
constexpr G4double PAD_X = 10.0*mm;
constexpr G4double PAD_Y = 10.0*mm;
constexpr G4double PAD_Z = 0.32*mm;
//Overall materials
//Thicknesses
constexpr G4double PAD_Absorber_Z = 3.5*mm;
constexpr G4double PAD_Glue1_Z  = 0.11*mm; //between W-Si
constexpr G4double PAD_Glue2_Z  = 0.13*mm; //between Si+FPC
constexpr G4double PAD_FPC_Z    = 0.28*mm;
constexpr G4double PAD_AirGap   = 1.0*mm;
constexpr G4double PAD_Layer_Thickness = PAD_Absorber_Z + PAD_Glue1_Z + PAD_Z + PAD_Glue2_Z + PAD_FPC_Z + PAD_AirGap;

//=================================================================
//                  Focal-H parameters
//=================================================================
//Dimentions of the detector
constexpr G4double HCAL_Z_Absorber = 30.0*mm;
constexpr G4double HCAL_Z_Scintillator = 2*mm;
constexpr G4double HCAL_Z_Gap = 0.0013*mm;
constexpr G4double HCal_Layer_Thickness = HCAL_Z_Absorber + HCAL_Z_Scintillator + HCAL_Z_Gap;

constexpr G4int HCALSiNumberOfLayers = 12;
constexpr G4double HCal_Si_Layer_Thickness = HCAL_Z_Absorber + PAD_Glue1_Z + PAD_Z + PAD_Glue2_Z + PAD_FPC_Z + PAD_AirGap;

constexpr G4int HCALSciNumberOfLayers = 30;
constexpr G4int HCALNumberOfTowersX = 6;
constexpr G4int HCALNumberOfTowersY = 6;       
constexpr G4int HCALNumberOfTowersZ = 2;
constexpr G4int NLayersHCALTower = HCALSciNumberOfLayers / HCALNumberOfTowersZ;

constexpr G4double HCAL_X_Tower = 100.0*mm; //25mm -> 10cm
constexpr G4double HCAL_Y_Tower = 100.0*mm;  
constexpr G4double HCAL_Tower_Gap = 20. *mm;    

#endif
