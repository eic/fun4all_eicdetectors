#include "PHG4ForwardDualReadoutDetector.h"
#include "PHG4ForwardDualReadoutDisplayAction.h"
#include "PHG4ForwardDualReadoutSteppingAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>              // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>               // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4PVParameterised.hh>
#include <Geant4/G4PVReplica.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4LogicalSkinSurface.hh>
#include <Geant4/G4LogicalBorderSurface.hh>

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________________
PHG4ForwardDualReadoutDetector::PHG4ForwardDualReadoutDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4ForwardDualReadoutDisplayAction*>(subsys->GetDisplayAction()))
  , m_SteppingAction(0)
  , _place_in_x(0.0 * mm)
  , _place_in_y(0.0 * mm)
  , _place_in_z(4000.0 * mm)
  , _center_offset_x(0.0 * mm)
  , _center_offset_y(0.0 * mm)
  , _quadratic_detector(0)
  , _rot_in_x(0.0)
  , _rot_in_y(0.0)
  , _rot_in_z(0.0)
  , _rMin1(50 * mm)
  , _rMax1(2620 * mm)
  , _rMin2(50 * mm)
  , _rMax2(3369 * mm)
  , _dZ(1000 * mm)
  , _sPhi(0)
  , _dPhi(2 * M_PI)
  , _tower_type(0)
  , _tower_readout(0.5 * mm)
  , _tower_dx(100 * mm)
  , _tower_dy(100 * mm)
  , _tower_dz(1000.0 * mm)
  , _scintFiber_diam(1.0 * mm)
  , _cerenkovFiber_diam(1.0 * mm)
  , _cerenkovFiber_material(0)
  , _tower_makeNotched(0)
  , _absorber_Material(0)
  , _materialScintillator("G4_POLYSTYRENE")
  , _materialAbsorber("G4_Fe")
  , _active(1)
  , _absorberactive(0)
  , _layer(0)
  , _blackhole(0)
  , _towerlogicnameprefix("hdrcaloTower")
  , _superdetector("NONE")
  , _mapping_tower_file("")
{
}
//_______________________________________________________________________
int PHG4ForwardDualReadoutDetector::IsInForwardDualReadout(G4VPhysicalVolume* volume) const
{
  if (volume->GetName().find(_towerlogicnameprefix) != string::npos)
  {
    if (volume->GetName().find("scintillator") != string::npos)
    {
      if (_active)
        return 1;
      else
        return 0;
    }
    else if (volume->GetName().find("cherenkov") != string::npos)
    {
      if (_active)
        return 1;
      else
        return 0;
    }
    //only record energy in actual absorber- drop energy lost in air gaps inside drcalo envelope
    else if (volume->GetName().find("absorber") != string::npos)
    {
      if (_absorberactive)
        return -1;
      else
        return 0;
    }
    else if (volume->GetName().find("envelope") != string::npos)
    {
      return 0;
    }
  }

  return 0;
}

//_______________________________________________________________________
void PHG4ForwardDualReadoutDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Begin Construction" << endl;
  }

  //Read parameters for detector construction from file
  ParseParametersFromTable();

  //Create the cone envelope = 'world volume' for the calorimeter
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* drcalo_envelope_solid;
  if(_quadratic_detector){
    // box with round cutout in the middle
    G4VSolid* beampipe_cutout = new G4Cons("hdrcalo_beampipe_cutout",
                                        0, _rMin1,
                                        0, _rMin1,
                                        _dZ / 2.0,
                                        _sPhi, _dPhi);
    drcalo_envelope_solid = new G4Box("hdrcalo_envelope_solid_precut",
                                        _rMax1,
                                        _rMax1,
                                        _tower_dz / 2.0);
    drcalo_envelope_solid = new G4SubtractionSolid(G4String("hdrcalo_envelope_solid"), drcalo_envelope_solid, beampipe_cutout
                                                            , 0 ,G4ThreeVector( 0 , 0 ,0.));
  } else {
    drcalo_envelope_solid = new G4Cons("hdrcalo_envelope_solid",
                                        _rMin1, _rMax1,
                                        _rMin2, _rMax2,
                                        _dZ / 2.0,
                                        _sPhi, _dPhi);
  }

  G4LogicalVolume* drcalo_envelope_log = new G4LogicalVolume(drcalo_envelope_solid, Air, G4String("hdrcalo_envelope"), 0, 0, 0);

  m_DisplayAction->AddVolume(drcalo_envelope_log, "FdrcaloEnvelope");

  //Define rotation attributes for envelope cone
  G4RotationMatrix drcalo_rotm;
  drcalo_rotm.rotateX(_rot_in_x);
  drcalo_rotm.rotateY(_rot_in_y);
  drcalo_rotm.rotateZ(_rot_in_z);

  //Place envelope cone in simulation
  ostringstream name_envelope;
  name_envelope.str("");
  name_envelope << _towerlogicnameprefix << "_envelope" << endl;

  new G4PVPlacement(G4Transform3D(drcalo_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z)),
                    drcalo_envelope_log, name_envelope.str().c_str(), logicWorld, 0, false, OverlapCheck());


  G4LogicalVolume* singletower;
  if(_tower_type==5) singletower = ConstructTowerFCStyle(0);
  else singletower = ConstructTower(0); //4x4 fibre tower with 2 scint and 2 Cherenkov fibres


  G4Material* material_air = G4Material::GetMaterial("G4_AIR");
 // number of towers in radial direction (on y axis)
  int rowNtow = (int) ( (_rMax1-(_tower_dy/2)) / _tower_dy);
  for(int row=rowNtow;row>=-rowNtow;row--){
    // pythagoras -> get available length in circular mother volume for towers
    // divide given length by tower width -> get number of towers that can be placed
    int currRowNtow = (int) ( ( 2* sqrt(pow(_rMax1,2)-pow( (abs(row)*_tower_dy) ,2)) ) / _tower_dy );
    if(currRowNtow==0) continue;
    // we want an odd number of towers to be symmetrically centered around 0
    if ( currRowNtow % 2 == 0) currRowNtow-=1;

    if( ( (abs(row)*_tower_dy) ) < _rMin1 ){ // _rMin1
      if(_center_offset_x!=0){
        // pythagoras -> get available length in circular mother volume for towers
        // divide given length by tower width -> get number of towers that can be placed
        int currRowNtowInner = (int) ( ( 2* sqrt(pow(_rMin1,2)-pow( (abs(row)*_tower_dy) ,2)) ) / _tower_dy );
        // we want an odd number of towers to be symmetrically centered around 0
        if(_quadratic_detector){
          if ( currRowNtowInner % 2 == 0) currRowNtowInner+=1;
          currRowNtowInner+=1;
        } else {
          if ( currRowNtowInner % 2 == 0) currRowNtowInner+=1;
        }
        int offsetrows = (int) ( _center_offset_x / _tower_dy );
        if ( offsetrows % 2 != 1) offsetrows-=1;

        int currRowNtowMod = currRowNtow;
        if(_quadratic_detector){
          currRowNtowMod = 2*rowNtow;
        }

        // create mother volume with space for currRowNtow towers along x-axis
        auto DRCalRowLeftSolid    = new G4Box("DRCalRowLeftBox" + std::to_string(row), ((currRowNtowMod - currRowNtowInner) / 2 + offsetrows) * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
        auto DRCalRowLeftLogical  = new G4LogicalVolume(DRCalRowLeftSolid,material_air,"DRCalRowLeftLogical" + std::to_string(row));
        // replicate singletower tower design currRowNtow times along x-axis
        new G4PVReplica("DRCalRowLeftPhysical" + std::to_string(row),singletower,DRCalRowLeftLogical,
                        kXAxis,((currRowNtowMod - currRowNtowInner) / 2 + offsetrows),_tower_dx);

        ostringstream name_row_twr_left;
        name_row_twr_left.str("");
        name_row_twr_left << _towerlogicnameprefix << "_row_" << row << "_left" << endl;
        new G4PVPlacement(0, G4ThreeVector( - ( ( currRowNtowInner / 2.0 ) * _tower_dx ) - ( ((currRowNtowMod - currRowNtowInner) / 2 - offsetrows) * _tower_dx / 2.0 ), (row*_tower_dy), 0),
                      DRCalRowLeftLogical, name_row_twr_left.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());

        // create mother volume with space for currRowNtow towers along x-axis
        auto DRCalRowRightSolid    = new G4Box("DRCalRowRightBox" + std::to_string(row), ((currRowNtowMod - currRowNtowInner) / 2 - offsetrows ) * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
        auto DRCalRowRightLogical  = new G4LogicalVolume(DRCalRowRightSolid,material_air,"DRCalRowRightLogical" + std::to_string(row));
        // replicate singletower tower design currRowNtow times along x-axis
        new G4PVReplica("DRCalRowRightPhysical" + std::to_string(row),singletower,DRCalRowRightLogical,
                        kXAxis,((currRowNtowMod - currRowNtowInner) / 2 - offsetrows ),_tower_dx);

        ostringstream name_row_twr_right;
        name_row_twr_right.str("");
        name_row_twr_right << _towerlogicnameprefix << "_row_" << row << "_right" << endl;
        new G4PVPlacement(0, G4ThreeVector( ( ( currRowNtowInner / 2.0 ) * _tower_dx ) + ( ((currRowNtowMod - currRowNtowInner) / 2 + offsetrows) * _tower_dx / 2.0 ), (row*_tower_dy), 0),
                      DRCalRowRightLogical, name_row_twr_right.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());
      } else {
        // pythagoras -> get available length in circular mother volume for towers
        // divide given length by tower width -> get number of towers that can be placed
        int currRowNtowInner = (int) ( ( 2* sqrt(pow(_rMin1,2)-pow( (abs(row)*_tower_dy) - (_tower_dy/2.0) ,2)) ) / _tower_dy );
        // we want an odd number of towers to be symmetrically centered around 0
        if ( currRowNtowInner % 2 != 0) currRowNtowInner-=1;
        // currRowNtowInner+=2;
        int currRowNtowMod = currRowNtow;
        if(_quadratic_detector){
          currRowNtowMod = rowNtow;
        }
        // create mother volume with space for currRowNtow towers along x-axis
        auto DRCalRowSolid    = new G4Box("DRCalRowBox" + std::to_string(row), (currRowNtowMod - currRowNtowInner) / 2 * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
        auto DRCalRowLogical  = new G4LogicalVolume(DRCalRowSolid,material_air,"DRCalRowLogical" + std::to_string(row));
        // replicate singletower tower design currRowNtow times along x-axis
        new G4PVReplica("DRCalRowPhysical" + std::to_string(row),singletower,DRCalRowLogical,
                        kXAxis,(currRowNtowMod - currRowNtowInner) / 2,_tower_dx);

        ostringstream name_row_twr;
        name_row_twr.str("");
        name_row_twr << _towerlogicnameprefix << "_row_" << row << "_left" << endl;
        new G4PVPlacement(0, G4ThreeVector( - ( ( currRowNtowInner / 2.0 ) * _tower_dx ) - ( (currRowNtowMod - currRowNtowInner) / 2 * _tower_dx / 2.0 ), (row*_tower_dy), 0),
                      DRCalRowLogical, name_row_twr.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());

        ostringstream name_row_twr2;
        name_row_twr2.str("");
        name_row_twr2 << _towerlogicnameprefix << "_row_" << row << "_left" << endl;
        new G4PVPlacement(0, G4ThreeVector( ( ( currRowNtowInner / 2.0 ) * _tower_dx ) + ( (currRowNtowMod - currRowNtowInner) / 2 * _tower_dx / 2.0 ), (row*_tower_dy), 0),
                      DRCalRowLogical, name_row_twr2.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());
      }

    } else {
      if(_quadratic_detector){
        // cout << currRowNtow << endl;
        // create mother volume with space for currRowNtow towers along x-axis
        auto DRCalRowSolid    = new G4Box("DRCalRowBox" + std::to_string(row), 2 * rowNtow * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
        auto DRCalRowLogical  = new G4LogicalVolume(DRCalRowSolid,material_air,"DRCalRowLogical" + std::to_string(row));
        // replicate singletower tower design currRowNtow times along x-axis
        new G4PVReplica("DRCalRowPhysical" + std::to_string(row),singletower,DRCalRowLogical,
                        kXAxis,2 * rowNtow,_tower_dx);

        ostringstream name_row_twr;
        name_row_twr.str("");
        name_row_twr << _towerlogicnameprefix << "_row_" << row << endl;
        new G4PVPlacement(0, G4ThreeVector(0, (row*_tower_dy), 0),
                      DRCalRowLogical, name_row_twr.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());
      } else {
        // create mother volume with space for currRowNtow towers along x-axis
        // cout << currRowNtow << endl;
        auto DRCalRowSolid    = new G4Box("DRCalRowBox" + std::to_string(row), currRowNtow * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
        auto DRCalRowLogical  = new G4LogicalVolume(DRCalRowSolid,material_air,"DRCalRowLogical" + std::to_string(row));
        // replicate singletower tower design currRowNtow times along x-axis
        new G4PVReplica("DRCalRowPhysical" + std::to_string(row),singletower,DRCalRowLogical,
                        kXAxis,currRowNtow,_tower_dx);

        ostringstream name_row_twr;
        name_row_twr.str("");
        name_row_twr << _towerlogicnameprefix << "_row_" << row << endl;
        new G4PVPlacement(0, G4ThreeVector(0, (row*_tower_dy), 0),
                      DRCalRowLogical, name_row_twr.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());
      }
    }
  }

  return;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4ForwardDualReadoutDetector::ConstructTower(int type)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Build logical volume for single tower..." << endl;
  }

  //create logical volume for single tower
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");
  // constructed tower base element
  G4VSolid* base_tower_solid = new G4Box(G4String("base_tower_solid"),
                                          _tower_dx / 2.0,
                                          _tower_dy / 2.0,
                                          _tower_dz / 2.0);

  G4LogicalVolume* base_tower_logic = new G4LogicalVolume(base_tower_solid,
                                                            material_air,
                                                            "base_tower_logic",
                                                            0, 0, 0);
  int maxsubtow = (int) ( (_tower_dx) / (_tower_readout));
  G4double addtowsize = (_tower_dx - (maxsubtow * _tower_readout))/maxsubtow;
  // 2x2 fiber tower base element
  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"),
                                          (_tower_readout + addtowsize) / 2.0,
                                          (_tower_readout + addtowsize) / 2.0,
                                          _tower_dz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  m_DisplayAction->AddVolume(single_tower_logic, "FdrcaloEnvelope");
  //create geometry volumes to place inside single_tower

  G4double diameter_fiber = _scintFiber_diam;
  G4double diameter_fiber_cherenkov = _cerenkovFiber_diam;
  G4double airgap = 0.1 * mm;

  // fibre cutout
  G4VSolid* single_cutout_tube  = new G4Tubs(G4String("single_cutout_tube"),
                                            0,
                                            ( diameter_fiber + airgap ) / 2.0,
                                            1.03 * _tower_dz / 1.0, //make it 1.03 times longer to ensure full cutout
                                            0.,2*M_PI*rad);
  G4VSolid* single_cutout_tube_cherenkov  = new G4Tubs(G4String("single_cutout_tube_cherenkov"),
                                            0,
                                            ( diameter_fiber_cherenkov + airgap ) / 2.0,
                                            1.03 * _tower_dz / 1.0, //make it 1.03 times longer to ensure full cutout
                                            0.,2*M_PI*rad);
  // notch cutout
  G4VSolid* single_cutout_box = new G4Box(G4String("single_cutout_box"),
                                          ( diameter_fiber + airgap ) / 2.0,
                                          1.03 * ( diameter_fiber + airgap ) / 4.0, //make it 1.03 times longer to ensure full cutout
                                          1.03 * _tower_dz / 1.0);
  G4VSolid* single_cutout_box_cherenkov = new G4Box(G4String("single_cutout_box_cherenkov"),
                                          ( diameter_fiber_cherenkov + airgap ) / 2.0,
                                          1.03 * ( diameter_fiber_cherenkov + airgap ) / 4.0, //make it 1.03 times longer to ensure full cutout
                                          1.03 * _tower_dz / 1.0);
  // absorber base object
  G4VSolid* solid_absorber_cher = new G4Box(G4String("solid_absorber_temp_cher"),
                                          (_tower_readout + addtowsize) / 4.0,
                                          (_tower_readout + addtowsize) / 4.0,
                                          _tower_dz / 2.0);
  G4VSolid* solid_absorber_scin = new G4Box(G4String("solid_absorber_temp_scin"),
                                          (_tower_readout + addtowsize) / 4.0,
                                          (_tower_readout + addtowsize) / 4.0,
                                          _tower_dz / 2.0);

  if(_tower_makeNotched){
    // cut out fiber hole
    solid_absorber_cher = new G4SubtractionSolid(G4String("solid_absorber_cher_f1"), solid_absorber_cher, single_cutout_tube_cherenkov
                                                              , 0 ,G4ThreeVector( 0 , ( (_tower_readout + addtowsize) / 4.0 ) - ( ( diameter_fiber_cherenkov + airgap ) / 2.0 ) ,0.)); // top right
    solid_absorber_scin = new G4SubtractionSolid(G4String("solid_absorber_scin_f1"), solid_absorber_scin, single_cutout_tube
                                                              , 0 ,G4ThreeVector( 0 , ( (_tower_readout + addtowsize) / 4.0 ) - ( ( diameter_fiber + airgap ) / 2.0 ) ,0.)); // top right
    // cut out notch
    solid_absorber_cher = new G4SubtractionSolid(G4String("solid_absorber_cher_box1"), solid_absorber_cher, single_cutout_box_cherenkov
                                                              , 0 ,G4ThreeVector( 0  , ( (_tower_readout + addtowsize) / 4.0 ) - ( ( diameter_fiber_cherenkov + airgap ) / 4.0 ) ,0.));
    solid_absorber_scin = new G4SubtractionSolid(G4String("solid_absorber_scin_box1"), solid_absorber_scin, single_cutout_box
                                                              , 0 ,G4ThreeVector( 0  , ( (_tower_readout + addtowsize) / 4.0 ) - ( ( diameter_fiber + airgap ) / 4.0 ) ,0.));
  } else {
    solid_absorber_cher = new G4SubtractionSolid(G4String("solid_absorber_temp_cher_f1"), solid_absorber_cher, single_cutout_tube_cherenkov
                                                              , 0 ,G4ThreeVector( 0 , 0 ,0.)); // top right
    solid_absorber_scin = new G4SubtractionSolid(G4String("solid_absorber_temp_scin_f1"), solid_absorber_scin, single_cutout_tube
                                                              , 0 ,G4ThreeVector( 0 , 0 ,0.)); // top right
  }
  G4VSolid* solid_scintillator  = new G4Tubs(G4String("single_scintillator_fiber"),
                                            0,
                                            diameter_fiber / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PI*rad);


  G4VSolid* solid_cherenkov  = new G4Tubs(G4String("single_cherenkov_fiber"),
                                            0,
                                            diameter_fiber_cherenkov / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PI*rad);


  G4Material* material_scintillator = GetScintillatorMaterial();


  G4NistManager* man = G4NistManager::Instance();
  G4Material* material_absorber;
  if(_absorber_Material==0)material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());
  else if(_absorber_Material==1)material_absorber = man->FindOrBuildMaterial("G4_W");
  else if(_absorber_Material==2)material_absorber = man->FindOrBuildMaterial("G4_Cu");
  else if(_absorber_Material==3)material_absorber = man->FindOrBuildMaterial("G4_Pb");
  else material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());


  G4LogicalVolume* logic_absorber_cher = new G4LogicalVolume(solid_absorber_cher,
                                                        material_absorber,
                                                        "absorber_solid_logic_cher",
                                                        0, 0, 0);

  G4LogicalVolume* logic_absorber_scin = new G4LogicalVolume(solid_absorber_scin,
                                                        material_absorber,
                                                        "absorber_solid_logic_scin",
                                                        0, 0, 0);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                    material_scintillator,
                                                    "hdrcalo_single_scintillator_fiber_logic",
                                                    0, 0, 0);

  G4Material *material_cherenkov;
  if(_cerenkovFiber_material==0) material_cherenkov = GetPMMAMaterial();
  else if(_cerenkovFiber_material==1) material_cherenkov = GetQuartzMaterial();
  else material_cherenkov = GetPMMAMaterial();

  G4LogicalVolume* logic_cherenk = new G4LogicalVolume(solid_cherenkov,
                                                    material_cherenkov,
                                                    "hdrcalo_single_cherenkov_fiber_logic",
                                                    0, 0, 0);

  m_DisplayAction->AddVolume(logic_absorber_cher, "Absorber");
  m_DisplayAction->AddVolume(logic_absorber_scin, "Absorber");

  m_DisplayAction->AddVolume(logic_scint, "Scintillator");
  m_DisplayAction->AddVolume(logic_cherenk, "Cherenkov");

  //place physical volumes for absorber and scintillator fiber

  ostringstream name_absorber;
  name_absorber.str("");
  name_absorber << _towerlogicnameprefix << "absorbersolid" << endl;

  ostringstream name_scintillator;
  name_scintillator.str("");
  name_scintillator << _towerlogicnameprefix << "singlescintillatorfiber" << endl;

  ostringstream name_cherenkov;
  name_cherenkov.str("");
  name_cherenkov << _towerlogicnameprefix << "_single_cherenkov_fiber"  << endl;


  new G4PVPlacement(0, G4ThreeVector( (_tower_readout + addtowsize) / 4.0,  (_tower_readout + addtowsize) / 4.0 , 0),
                    logic_absorber_cher,
                    name_absorber.str().c_str()+std::to_string(1),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  new G4PVPlacement(0, G4ThreeVector( (_tower_readout + addtowsize) / 4.0,  -(_tower_readout + addtowsize) / 4.0 , 0),
                    logic_absorber_scin,
                    name_absorber.str().c_str()+std::to_string(2),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  new G4PVPlacement(0, G4ThreeVector( -(_tower_readout + addtowsize) / 4.0,  (_tower_readout + addtowsize) / 4.0 , 0),
                    logic_absorber_scin,
                    name_absorber.str().c_str()+std::to_string(3),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  new G4PVPlacement(0, G4ThreeVector( -(_tower_readout + addtowsize) / 4.0,  -(_tower_readout + addtowsize) / 4.0 , 0),
                    logic_absorber_cher,
                    name_absorber.str().c_str()+std::to_string(4),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  if(_tower_makeNotched){
    // place scintillator fibers (top left, bottom right)
    new G4PVPlacement(0, G4ThreeVector( -(_tower_readout + addtowsize) / 4.0,  ( (_tower_readout + addtowsize) / 2.0 ) - ( ( diameter_fiber + airgap ) / 2.0 ) , 0),
                      logic_scint,
                      name_scintillator.str().c_str()+std::to_string(1),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( (_tower_readout + addtowsize) / 4.0,  - ( ( diameter_fiber + airgap ) / 2.0 ) , 0),
                      logic_scint,
                      name_scintillator.str().c_str()+std::to_string(2),
                      single_tower_logic,
                      0, 0, OverlapCheck());

    // place cherenkov fibers (top right, bottom left)
    new G4PVPlacement(0, G4ThreeVector( (_tower_readout + addtowsize) / 4.0 ,  ( (_tower_readout + addtowsize) / 2.0 ) - ( ( diameter_fiber_cherenkov + airgap ) / 2.0 ) , 0),
                      logic_cherenk,
                      name_cherenkov.str().c_str()+std::to_string(1),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( -(_tower_readout + addtowsize) / 4.0 ,  - ( ( diameter_fiber_cherenkov + airgap ) / 2.0 ) , 0),
                      logic_cherenk,
                      name_cherenkov.str().c_str()+std::to_string(2),
                      single_tower_logic,
                      0, 0, OverlapCheck());
  } else {
    // place scintillator fibers (top left, bottom right)
    new G4PVPlacement(0, G4ThreeVector( -(_tower_readout + addtowsize) / 4.0,  (_tower_readout + addtowsize) / 4.0 , 0),
                      logic_scint,
                      name_scintillator.str().c_str()+std::to_string(1),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( (_tower_readout + addtowsize) / 4.0,  - (_tower_readout + addtowsize) / 4.0 , 0),
                      logic_scint,
                      name_scintillator.str().c_str()+std::to_string(2),
                      single_tower_logic,
                      0, 0, OverlapCheck());

    // place cherenkov fibers (top right, bottom left)
    new G4PVPlacement(0, G4ThreeVector( (_tower_readout + addtowsize) / 4.0 ,  (_tower_readout + addtowsize) / 4.0 , 0),
                      logic_cherenk,
                      name_cherenkov.str().c_str()+std::to_string(3),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( -(_tower_readout + addtowsize) / 4.0 ,  - (_tower_readout + addtowsize) / 4.0 , 0),
                      logic_cherenk,
                      name_cherenkov.str().c_str()+std::to_string(4),
                      single_tower_logic,
                      0, 0, OverlapCheck());
  }


  int rowNtow = (int) ( (_tower_dx) / (_tower_readout + addtowsize));
  for(int row=(rowNtow / 2);row>=-rowNtow / 2;row--){
      // create mother volume with space for currRowNtow towers along x-axis
      auto DRCalRowSolid    = new G4Box("DRCalRowBoxBase" + std::to_string(row), _tower_dx / 2.0, (_tower_readout + addtowsize) / 2.0, _tower_dz / 2.0);
      auto DRCalRowLogical  = new G4LogicalVolume(DRCalRowSolid,material_air,"DRCalRowLogicalBase" + std::to_string(row));
      // replicate singletower tower design currRowNtow times along x-axis
      new G4PVReplica("DRCalRowPhysicalBase" + std::to_string(row),single_tower_logic,DRCalRowLogical,
                      kXAxis, rowNtow, _tower_readout + addtowsize);

      m_DisplayAction->AddVolume(DRCalRowLogical, "FdrcaloEnvelope");
      ostringstream name_row_twr;
      name_row_twr.str("");
      name_row_twr << _towerlogicnameprefix << "_row_" << row << endl;
      new G4PVPlacement(0, G4ThreeVector(0, (row * (_tower_readout + addtowsize)), 0),
                    DRCalRowLogical, name_row_twr.str().c_str(), base_tower_logic, 0, false, OverlapCheck());
  }

  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Building logical volume for single tower done." << endl;
  }

  return base_tower_logic;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4ForwardDualReadoutDetector::ConstructTowerFCStyle(int type)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Build logical volume for single tower..." << endl;
  }

  //create logical volume for single tower
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");
  // 2x2 tower base element
  G4VSolid* base_tower_solid = new G4Box(G4String("base_tower_solid"),
                                          _tower_dx / 2.0,
                                          _tower_dy / 2.0,
                                          _tower_dz / 2.0);

  G4LogicalVolume* base_tower_logic = new G4LogicalVolume(base_tower_solid,
                                                            material_air,
                                                            "base_tower_logic",
                                                            0, 0, 0);
  G4double copperTubeDiam = _tower_readout / 2;
  int maxsubtow = (int) ( (_tower_dx) / (2 * copperTubeDiam));
  G4double addtowsize = (_tower_dx - (maxsubtow * 2 * copperTubeDiam))/maxsubtow;
  // 2x2 tower base element
  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"),
                                          (2 * copperTubeDiam + addtowsize) / 2.0,
                                          (2 * copperTubeDiam + addtowsize) / 2.0,
                                          _tower_dz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  m_DisplayAction->AddVolume(single_tower_logic, "FdrcaloEnvelope");
  //create geometry volumes to place inside single_tower

  G4double diameter_fiber = _scintFiber_diam;
  G4double diameter_fiber_cherenkov = _cerenkovFiber_diam;

  // fibre cutout
  G4VSolid* solid_absorber  = new G4Tubs(G4String("ttl_copper_tube_solid"),
                                            0,
                                            copperTubeDiam / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PI*rad);

  G4VSolid* solid_scintillator  = new G4Tubs(G4String("single_scintillator_fiber"),
                                            0,
                                            diameter_fiber / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PI*rad);


  G4VSolid* solid_cherenkov  = new G4Tubs(G4String("single_cherenkov_fiber"),
                                            0,
                                            diameter_fiber_cherenkov / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PI*rad);


  G4Material* material_scintillator = GetScintillatorMaterial();


  G4NistManager* man = G4NistManager::Instance();
  G4Material* material_absorber;
  if(_absorber_Material==0)material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());
  else if(_absorber_Material==1)material_absorber = man->FindOrBuildMaterial("G4_W");
  else if(_absorber_Material==2)material_absorber = man->FindOrBuildMaterial("G4_Cu");
  else if(_absorber_Material==3)material_absorber = man->FindOrBuildMaterial("G4_Pb");
  else material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());


  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "absorber_solid_logic",
                                                        0, 0, 0);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                    material_scintillator,
                                                    "hdrcalo_single_scintillator_fiber_logic",
                                                    0, 0, 0);

  G4Material *material_cherenkov;
  if(_cerenkovFiber_material==0) material_cherenkov = GetPMMAMaterial();
  else if(_cerenkovFiber_material==1) material_cherenkov = GetQuartzMaterial();
  else material_cherenkov = GetPMMAMaterial();

  G4LogicalVolume* logic_cherenk = new G4LogicalVolume(solid_cherenkov,
                                                    material_cherenkov,
                                                    "hdrcalo_single_cherenkov_fiber_logic",
                                                    0, 0, 0);

  m_DisplayAction->AddVolume(logic_absorber, "Absorber");
  m_DisplayAction->AddVolume(logic_scint, "Scintillator");
  m_DisplayAction->AddVolume(logic_cherenk, "Cherenkov");

  //place physical volumes for absorber and scintillator fiber

  ostringstream name_absorber;
  name_absorber.str("");
  name_absorber << _towerlogicnameprefix << "absorbersolid" << endl;

  ostringstream name_scintillator;
  name_scintillator.str("");
  name_scintillator << _towerlogicnameprefix << "singlescintillatorfiber" << endl;

  ostringstream name_cherenkov;
  name_cherenkov.str("");
  name_cherenkov << _towerlogicnameprefix << "_single_cherenkov_fiber"  << endl;

  // place copper rods
  new G4PVPlacement(0, G4ThreeVector( -copperTubeDiam/2,  copperTubeDiam/2 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( -copperTubeDiam/2,  -copperTubeDiam/2 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( copperTubeDiam/2,  copperTubeDiam/2 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( copperTubeDiam/2,  -copperTubeDiam/2 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  // place scintillator fibers (top left, bottom right)
  new G4PVPlacement(0, G4ThreeVector( 0,  0 , 0),
                    logic_scint,
                    name_scintillator.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( copperTubeDiam, copperTubeDiam , 0),
                    logic_scint,
                    name_scintillator.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  // place cherenkov fibers (top right, bottom left)
  new G4PVPlacement(0, G4ThreeVector( copperTubeDiam, 0 , 0),
                    logic_cherenk,
                    name_cherenkov.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( 0, copperTubeDiam , 0),
                    logic_cherenk,
                    name_cherenkov.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());


  int rowNtow = (int) ( (_tower_dx) / (2 * copperTubeDiam + addtowsize));
  for(int row=(rowNtow / 2);row>=-rowNtow / 2;row--){
      // create mother volume with space for currRowNtow towers along x-axis
      auto DRCalRowSolid    = new G4Box("DRCalRowBox", _tower_dx / 2.0, (2 * copperTubeDiam + addtowsize) / 2.0, _tower_dz / 2.0);
      auto DRCalRowLogical  = new G4LogicalVolume(DRCalRowSolid,material_air,"DRCalRowLogical");
      // replicate singletower tower design currRowNtow times along x-axis
      new G4PVReplica("DRCalRowPhysical",single_tower_logic,DRCalRowLogical,
                      kXAxis, rowNtow, 2 * copperTubeDiam + addtowsize);

      m_DisplayAction->AddVolume(DRCalRowLogical, "FdrcaloEnvelope");
      ostringstream name_row_twr;
      name_row_twr.str("");
      name_row_twr << _towerlogicnameprefix << "_row_" << row << endl;
      new G4PVPlacement(0, G4ThreeVector(0, (row * (2 * copperTubeDiam + addtowsize)), 0),
                    DRCalRowLogical, name_row_twr.str().c_str(), base_tower_logic, 0, false, OverlapCheck());
  }

  // new G4PVPlacement(0, G4ThreeVector( 0, 0 , 0),
  //                   single_tower_logic,
  //                   name_cherenkov.str().c_str(),
  //                   base_tower_logic,
  //                   0, 0, OverlapCheck());

  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Building logical volume for single tower done." << endl;
  }

  return base_tower_logic;
}

//_______________________________________________________________________
G4Material*
PHG4ForwardDualReadoutDetector::GetScintillatorMaterial()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Making Scintillator material..." << endl;
  }


  G4MaterialPropertiesTable *tab = new G4MaterialPropertiesTable();

  const G4int ntab = 31;
  tab->AddConstProperty("FASTTIMECONSTANT", 2.8*ns); // was 6
  // tab->AddConstProperty("SCINTILLATIONYIELD", 13.9/keV); // was 200/MEV nominal  10
  tab->AddConstProperty("SCINTILLATIONYIELD", 200/MeV); // was 200/MEV nominal, should maybe be 13.9/keV
  tab->AddConstProperty("RESOLUTIONSCALE", 1.0);

  G4double opt_en[] =
    { 1.37760*eV, 1.45864*eV, 1.54980*eV, 1.65312*eV, 1.71013*eV, 1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV,
      2.10143*eV, 2.13766*eV, 2.17516*eV, 2.21400*eV, 2.25426*eV, 2.29600*eV, 2.33932*eV, 2.38431*eV, 2.43106*eV, 2.47968*eV,
      2.53029*eV, 2.58300*eV, 2.63796*eV, 2.69531*eV, 2.75520*eV, 2.81782*eV, 2.88335*eV, 2.95200*eV, 3.09960*eV, 3.54241*eV,
      4.13281*eV }; // 350 - 800 nm
  G4double scin_fast[] =
    { 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0003, 0.0008, 0.0032,
      0.0057, 0.0084, 0.0153, 0.0234, 0.0343, 0.0604, 0.0927, 0.1398, 0.2105, 0.2903,
      0.4122, 0.5518, 0.7086, 0.8678, 1.0000, 0.8676, 0.2311, 0.0033, 0.0012, 0.0000,
      0 };
  tab->AddProperty("FASTCOMPONENT", opt_en, scin_fast, ntab);

  G4double opt_r[] =
    { 1.5749, 1.5764, 1.5782, 1.5803, 1.5815, 1.5829, 1.5845, 1.5862, 1.5882, 1.5904,
      1.5914, 1.5924, 1.5935, 1.5947, 1.5959, 1.5972, 1.5986, 1.6000, 1.6016, 1.6033,
      1.6051, 1.6070, 1.6090, 1.6112, 1.6136, 1.6161, 1.6170, 1.6230, 1.62858, 1.65191,
      1.69165 };
  tab->AddProperty("RINDEX", opt_en, opt_r, ntab);

  G4double opt_abs[] =
    { 2.714*m, 3.619*m, 5.791*m, 4.343*m, 7.896*m, 5.429*m, 36.19*m, 17.37*m, 36.19*m, 5.429*m,
      13.00*m, 14.50*m, 16.00*m, 18.00*m, 16.50*m, 17.00*m, 14.00*m, 16.00*m, 15.00*m, 14.50*m,
      13.00*m, 12.00*m, 10.00*m, 8.000*m, 7.238*m, 4.000*m, 1.200*m, 0.500*m, 0.200*m, 0.200*m,
      0.100*m };
  tab->AddProperty("ABSLENGTH", opt_en, opt_abs, ntab);

	G4double density;
	G4int ncomponents;
  // G4Material* material_G4_POLYSTYRENE = G4Material::GetMaterial(_materialScintillator.c_str());
  G4Material* material_G4_POLYSTYRENE = new G4Material("G4_POLYSTYRENE", density = 1.05 * g / cm3, ncomponents = 2);
  material_G4_POLYSTYRENE->AddElement(G4Element::GetElement("C"), 8);
  material_G4_POLYSTYRENE->AddElement(G4Element::GetElement("H"), 8);
  material_G4_POLYSTYRENE->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  material_G4_POLYSTYRENE->SetMaterialPropertiesTable(tab);

  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector:  Making Scintillator material done." << endl;
  }

  return material_G4_POLYSTYRENE;
}

//_______________________________________________________________________
G4Material*
PHG4ForwardDualReadoutDetector::GetPMMAMaterial()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Making PMMA material..." << endl;
  }

	G4double density;
	G4int ncomponents;

  G4Material* material_PMMA = new G4Material("PMMA", density = 1.18 * g / cm3, ncomponents = 3);
  material_PMMA->AddElement(G4Element::GetElement("C"), 5);
  material_PMMA->AddElement(G4Element::GetElement("H"), 8);
  material_PMMA->AddElement(G4Element::GetElement("O"), 2);

  const G4int nEntries = 31;

  G4double photonEnergy[nEntries] =
    { 1.37760*eV, 1.45864*eV, 1.54980*eV, 1.65312*eV, 1.71013*eV, 1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV,
      2.10143*eV, 2.13766*eV, 2.17516*eV, 2.21400*eV, 2.25426*eV, 2.29600*eV, 2.33932*eV, 2.38431*eV, 2.43106*eV, 2.47968*eV,
      2.53029*eV, 2.58300*eV, 2.63796*eV, 2.69531*eV, 2.75520*eV, 2.81782*eV, 2.88335*eV, 2.95200*eV, 3.09960*eV, 3.54241*eV,
      4.13281*eV };
  G4double refractiveIndexWLSfiber[nEntries] =
    { 1.4852, 1.4859, 1.4867, 1.4877, 1.4882, 1.4888, 1.4895, 1.4903, 1.4911, 1.4920,
      1.4924, 1.4929, 1.4933, 1.4938, 1.4943, 1.4948, 1.4954, 1.4960, 1.4966, 1.4973,
      1.4981, 1.4989, 1.4997, 1.5006, 1.5016, 1.5026, 1.5038, 1.5050, 1.5052, 1.5152,
      1.5306 };

  G4double absWLSfiber[nEntries] =
    { 0.414*m, 0.965*m, 2.171*m, 4.343*m, 1.448*m, 4.343*m, 14.48*m, 21.71*m, 8.686*m, 39.48*m,
      48.25*m, 54.29*m, 57.91*m, 54.29*m, 33.40*m, 31.02*m, 43.43*m, 43.43*m, 41.36*m, 39.48*m,
      37.76*m, 36.19*m, 36.19*m, 33.40*m, 31.02*m, 28.95*m, 25.55*m, 24.13*m, 21.71*m, 2.171*m,
      0.434*m };


  // Add entries into properties table
  G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
  mptWLSfiber->AddProperty("RINDEX",photonEnergy,refractiveIndexWLSfiber,nEntries);
  mptWLSfiber->AddProperty("ABSLENGTH",photonEnergy,absWLSfiber,nEntries);
  material_PMMA->SetMaterialPropertiesTable(mptWLSfiber);
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector:  Making PMMA material done." << endl;
  }

  return material_PMMA;
}

//_______________________________________________________________________
G4Material*
PHG4ForwardDualReadoutDetector::GetQuartzMaterial()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Making Scintillator material..." << endl;
  }

  G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();

	G4double density;
	G4int ncomponents;

  // G4Material* material_Quartz = G4Material::GetMaterial("Quartz");
  G4Material *material_Quartz = new G4Material("Quartz", density = 2.200 * g / cm3, ncomponents = 2);
  material_Quartz->AddElement(G4Element::GetElement("Si"), 1);
  material_Quartz->AddElement(G4Element::GetElement("O"), 2);

  const G4int nEntriesQuartz = 279;

  G4double photonEnergyQuartz[nEntriesQuartz] =
    { 5.9040*eV, 5.7667*eV, 5.6356*eV, 5.5104*eV, 5.3906*eV, 5.2759*eV, 5.1660*eV, 5.0606*eV, 4.9594*eV, 4.8621*eV,
      4.7686*eV, 4.6786*eV, 4.5920*eV, 4.5085*eV, 4.4280*eV, 4.3503*eV, 4.2753*eV, 4.2029*eV, 4.1328*eV, 4.0651*eV,
      3.9995*eV, 3.9360*eV, 3.8745*eV, 3.8149*eV, 3.7571*eV, 3.7010*eV, 3.6466*eV, 3.5937*eV, 3.5424*eV, 3.4925*eV,
      3.4440*eV, 3.3968*eV, 3.3509*eV, 3.3062*eV, 3.2627*eV, 3.2204*eV, 3.1791*eV, 3.1388*eV, 3.0996*eV, 3.0613*eV,
      3.0240*eV, 2.9876*eV, 2.9520*eV, 2.9173*eV, 2.8834*eV, 2.8502*eV, 2.8178*eV, 2.7862*eV, 2.7552*eV, 2.7249*eV,
      2.6953*eV, 2.6663*eV, 2.6380*eV, 2.6102*eV, 2.5830*eV, 2.5564*eV, 2.5303*eV, 2.5047*eV, 2.4797*eV, 2.4551*eV,
      2.4311*eV, 2.4075*eV, 2.3843*eV, 2.3616*eV, 2.3393*eV, 2.3175*eV, 2.2960*eV, 2.2749*eV, 2.2543*eV, 2.2339*eV,
      2.2140*eV, 2.1944*eV, 2.1752*eV, 2.1562*eV, 2.1377*eV, 2.1194*eV, 2.1014*eV, 2.0838*eV, 2.0664*eV, 2.0493*eV,
      2.0325*eV, 2.0160*eV, 1.9997*eV, 1.9837*eV, 1.9680*eV, 1.9525*eV, 1.9373*eV, 1.9222*eV, 1.9074*eV, 1.8929*eV,
      1.8785*eV, 1.8644*eV, 1.8505*eV, 1.8368*eV, 1.8233*eV, 1.8100*eV, 1.7969*eV, 1.7839*eV, 1.7712*eV, 1.7463*eV,
      1.7220*eV, 1.6984*eV, 1.6755*eV, 1.6531*eV, 1.6314*eV, 1.6102*eV, 1.5895*eV, 1.5694*eV, 1.5498*eV, 1.5307*eV,
      1.5120*eV, 1.4938*eV, 1.4760*eV, 1.4586*eV, 1.4417*eV, 1.4251*eV, 1.4089*eV, 1.3931*eV, 1.3776*eV, 1.3625*eV,
      1.3477*eV, 1.3332*eV, 1.3190*eV, 1.3051*eV, 1.2915*eV, 1.2782*eV, 1.2651*eV, 1.2524*eV, 1.2398*eV, 1.2276*eV,
      1.2155*eV, 1.2037*eV, 1.1922*eV, 1.1808*eV, 1.1697*eV, 1.1587*eV, 1.1480*eV, 1.1375*eV, 1.1271*eV, 1.1170*eV,
      1.1070*eV, 1.0972*eV, 1.0876*eV, 1.0781*eV, 1.0688*eV, 1.0597*eV, 1.0507*eV, 1.0419*eV, 1.0332*eV, 1.0247*eV,
      1.0163*eV, 1.0080*eV, 0.9999*eV, 0.9919*eV, 0.9840*eV, 0.9763*eV, 0.9686*eV, 0.9611*eV, 0.9537*eV, 0.9464*eV,
      0.9393*eV, 0.9322*eV, 0.9253*eV, 0.9184*eV, 0.9116*eV, 0.9050*eV, 0.8984*eV, 0.8920*eV, 0.8856*eV, 0.8793*eV,
      0.8731*eV, 0.8670*eV, 0.8610*eV, 0.8551*eV, 0.8492*eV, 0.8434*eV, 0.8377*eV, 0.8321*eV, 0.8266*eV, 0.8211*eV,
      0.8157*eV, 0.8104*eV, 0.8051*eV, 0.7999*eV, 0.7948*eV, 0.7897*eV, 0.7847*eV, 0.7798*eV, 0.7749*eV, 0.7701*eV,
      0.7653*eV, 0.7606*eV, 0.7560*eV, 0.7514*eV, 0.7469*eV, 0.7424*eV, 0.7380*eV, 0.7336*eV, 0.7293*eV, 0.7251*eV,
      0.7208*eV, 0.7167*eV, 0.7126*eV, 0.7085*eV, 0.7045*eV, 0.7005*eV, 0.6965*eV, 0.6926*eV, 0.6888*eV, 0.6850*eV,
      0.6812*eV, 0.6775*eV, 0.6738*eV, 0.6702*eV, 0.6666*eV, 0.6630*eV, 0.6595*eV, 0.6560*eV, 0.6525*eV, 0.6491*eV,
      0.6458*eV, 0.6424*eV, 0.6391*eV, 0.6358*eV, 0.6326*eV, 0.6294*eV, 0.6262*eV, 0.6230*eV, 0.6199*eV, 0.6168*eV,
      0.6138*eV, 0.6108*eV, 0.6078*eV, 0.6048*eV, 0.6019*eV, 0.5990*eV, 0.5961*eV, 0.5932*eV, 0.5904*eV, 0.5876*eV,
      0.5848*eV, 0.5821*eV, 0.5794*eV, 0.5767*eV, 0.5740*eV, 0.5714*eV, 0.5687*eV, 0.5661*eV, 0.5636*eV, 0.5610*eV,
      0.5585*eV, 0.5560*eV, 0.5535*eV, 0.5510*eV, 0.5486*eV, 0.5462*eV, 0.5438*eV, 0.5414*eV, 0.5391*eV, 0.5367*eV,
      0.5344*eV, 0.5321*eV, 0.5298*eV, 0.5276*eV, 0.5254*eV, 0.5231*eV, 0.5209*eV, 0.5188*eV, 0.5166*eV, 0.5145*eV,
      0.5123*eV, 0.5102*eV, 0.5081*eV, 0.5061*eV, 0.5040*eV, 0.5020*eV, 0.4999*eV, 0.4979*eV, 0.4959*eV};
  G4double refractiveIndexQuartz[nEntriesQuartz] =
    { 1.5384, 1.5332, 1.5285, 1.5242, 1.5202, 1.5166, 1.5133, 1.5103, 1.5074, 1.5048,
      1.5024, 1.5001, 1.4980, 1.4960, 1.4942, 1.4924, 1.4908, 1.4892, 1.4878, 1.4864,
      1.4851, 1.4839, 1.4827, 1.4816, 1.4806, 1.4796, 1.4787, 1.4778, 1.4769, 1.4761,
      1.4753, 1.4745, 1.4738, 1.4731, 1.4725, 1.4719, 1.4713, 1.4707, 1.4701, 1.4696,
      1.4691, 1.4686, 1.4681, 1.4676, 1.4672, 1.4668, 1.4663, 1.4660, 1.4656, 1.4652,
      1.4648, 1.4645, 1.4641, 1.4638, 1.4635, 1.4632, 1.4629, 1.4626, 1.4623, 1.4621,
      1.4618, 1.4615, 1.4613, 1.4610, 1.4608, 1.4606, 1.4603, 1.4601, 1.4599, 1.4597,
      1.4595, 1.4593, 1.4591, 1.4589, 1.4587, 1.4586, 1.4584, 1.4582, 1.4580, 1.4579,
      1.4577, 1.4576, 1.4574, 1.4572, 1.4571, 1.4570, 1.4568, 1.4567, 1.4565, 1.4564,
      1.4563, 1.4561, 1.4560, 1.4559, 1.4558, 1.4556, 1.4555, 1.4554, 1.4553, 1.4551,
      1.4549, 1.4546, 1.4544, 1.4542, 1.4540, 1.4539, 1.4537, 1.4535, 1.4533, 1.4531,
      1.4530, 1.4528, 1.4527, 1.4525, 1.4523, 1.4522, 1.4520, 1.4519, 1.4518, 1.4516,
      1.4515, 1.4513, 1.4512, 1.4511, 1.4509, 1.4508, 1.4507, 1.4505, 1.4504, 1.4503,
      1.4502, 1.4500, 1.4499, 1.4498, 1.4497, 1.4496, 1.4494, 1.4493, 1.4492, 1.4491,
      1.4490, 1.4489, 1.4487, 1.4486, 1.4485, 1.4484, 1.4483, 1.4482, 1.4481, 1.4479,
      1.4478, 1.4477, 1.4476, 1.4475, 1.4474, 1.4473, 1.4471, 1.4470, 1.4469, 1.4468,
      1.4467, 1.4466, 1.4465, 1.4464, 1.4462, 1.4461, 1.4460, 1.4459, 1.4458, 1.4457,
      1.4455, 1.4454, 1.4453, 1.4452, 1.4451, 1.4450, 1.4449, 1.4447, 1.4446, 1.4445,
      1.4444, 1.4443, 1.4441, 1.4440, 1.4439, 1.4438, 1.4437, 1.4435, 1.4434, 1.4433,
      1.4432, 1.4431, 1.4429, 1.4428, 1.4427, 1.4426, 1.4424, 1.4423, 1.4422, 1.4420,
      1.4419, 1.4418, 1.4417, 1.4415, 1.4414, 1.4413, 1.4411, 1.4410, 1.4409, 1.4407,
      1.4406, 1.4405, 1.4403, 1.4402, 1.4401, 1.4399, 1.4398, 1.4397, 1.4395, 1.4394,
      1.4392, 1.4391, 1.4389, 1.4388, 1.4387, 1.4385, 1.4384, 1.4382, 1.4381, 1.4379,
      1.4378, 1.4376, 1.4375, 1.4373, 1.4372, 1.4370, 1.4369, 1.4367, 1.4366, 1.4364,
      1.4363, 1.4361, 1.4360, 1.4358, 1.4357, 1.4355, 1.4353, 1.4352, 1.4350, 1.4349,
      1.4347, 1.4345, 1.4344, 1.4342, 1.4340, 1.4339, 1.4337, 1.4335, 1.4334, 1.4332,
      1.4330, 1.4328, 1.4327, 1.4325, 1.4323, 1.4322, 1.4320, 1.4318, 1.4316, 1.4314,
      1.4313, 1.4311, 1.4309, 1.4307, 1.4305, 1.4304, 1.4302, 1.4300, 1.4298 };
  mptWLSfiber->AddProperty("RINDEX",photonEnergyQuartz,refractiveIndexQuartz,nEntriesQuartz);


  const G4int nEntries_Quartz = 14;
  G4double PhotonEnergy_Quartz[nEntries_Quartz] =
    { 2.21*eV, 2.30*eV, 2.38*eV, 2.48*eV, 2.58*eV, 2.70*eV, 2.82*eV, 2.95*eV, 3.10*eV, 3.26*eV,
      3.44*eV, 3.65*eV, 3.88*eV, 4.13*eV };
  G4double Quartz_Abs[nEntries_Quartz] =
    { 550.7*mm, 530.7*mm, 590.1*mm, 490.7*mm, 470.7*mm, 520.3*mm, 500.0*mm, 470.7*mm, 450.5*mm, 270.5*mm,
      190.1*mm,  60.9*mm,  10.6*mm,   4.0*mm};
  mptWLSfiber->AddProperty("ABSLENGTH",  PhotonEnergy_Quartz, Quartz_Abs,  nEntries_Quartz);

  material_Quartz->SetMaterialPropertiesTable(mptWLSfiber);

  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector:  Making Scintillator material done." << endl;
  }

  return material_Quartz;
}


int PHG4ForwardDualReadoutDetector::ParseParametersFromTable()
{
  //Open the datafile, if it won't open return an error
  ifstream istream_mapping;
  istream_mapping.open(_mapping_tower_file);
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in PHG4ForwardHcalDetector: Failed to open mapping file " << _mapping_tower_file << std::endl;
    gSystem->Exit(1);
  }

  //loop over lines in file
  string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {
    //Skip lines starting with / including a '#'
    if (line_mapping.find("#") != string::npos)
    {
      if (Verbosity() > 0)
      {
        std::cout << "PHG4ForwardHcalDetector: SKIPPING line in mapping file: " << line_mapping << std::endl;
      }
      continue;
    }

    istringstream iss(line_mapping);
      //If this line is not a comment and not a tower, save parameter as string / value.
      string parname;
      G4double parval;

      //read string- break if error
      if (!(iss >> parname >> parval))
      {
        cout << "ERROR in PHG4ForwardHcalDetector: Failed to read line in mapping file " << _mapping_tower_file << std::endl;
        gSystem->Exit(1);
      }

      m_GlobalParameterMap.insert(make_pair(parname, parval));
  }

  //Update member variables for global parameters based on parsed parameter file
  std::map<string, G4double>::iterator parit;


  parit = m_GlobalParameterMap.find("Gtype");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_type = parit->second;
  }

  parit = m_GlobalParameterMap.find("Gtower_readout");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_readout = parit->second * cm;
    m_SteppingAction->SetTowerReadout(_tower_readout);
  }

  parit = m_GlobalParameterMap.find("Gtower_dx");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_dx = parit->second * cm;
    m_SteppingAction->SetTowerSize(_tower_dx);
  }

  parit = m_GlobalParameterMap.find("Gtower_dy");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_dy = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gtower_dz");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_dz = parit->second * cm;
  }

  // new start
  parit = m_GlobalParameterMap.find("Scint_Diam");
  if (parit != m_GlobalParameterMap.end())
  {
    _scintFiber_diam = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Cerenkov_Diam");
  if (parit != m_GlobalParameterMap.end())
  {
    _cerenkovFiber_diam = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Cerenkov_Material");
  if (parit != m_GlobalParameterMap.end())
  {
    _cerenkovFiber_material = parit->second;
  }

  parit = m_GlobalParameterMap.find("NotchCutout");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_makeNotched = parit->second;
  }

  parit = m_GlobalParameterMap.find("Absorber_Material");
  if (parit != m_GlobalParameterMap.end())
  {
    _absorber_Material = parit->second;
  }
  // new end

  parit = m_GlobalParameterMap.find("Gr1_inner");
  if (parit != m_GlobalParameterMap.end())
  {
    _rMin1 = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gr1_outer");
  if (parit != m_GlobalParameterMap.end())
  {
    _rMax1 = parit->second * cm;
    m_SteppingAction->SetDetectorSize(_rMax1);
  }

  parit = m_GlobalParameterMap.find("Gr2_inner");
  if (parit != m_GlobalParameterMap.end())
  {
    _rMin2 = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gr2_outer");
  if (parit != m_GlobalParameterMap.end())
  {
    _rMax2 = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gdz");
  if (parit != m_GlobalParameterMap.end())
  {
    _dZ = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gx0");
  if (parit != m_GlobalParameterMap.end())
  {
    _place_in_x = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gy0");
  if (parit != m_GlobalParameterMap.end())
  {
    _place_in_y = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gz0");
  if (parit != m_GlobalParameterMap.end())
  {
    _place_in_z = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Center_Offset_x");
  if (parit != m_GlobalParameterMap.end())
  {
    _center_offset_x = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Center_Offset_y");
  if (parit != m_GlobalParameterMap.end())
  {
    _center_offset_y = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Quadratic_Detector");
  if (parit != m_GlobalParameterMap.end())
  {
    _quadratic_detector = parit->second;
  }

  parit = m_GlobalParameterMap.find("Grot_x");
  if (parit != m_GlobalParameterMap.end())
  {
    _rot_in_x = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Grot_y");
  if (parit != m_GlobalParameterMap.end())
  {
    _rot_in_y = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Grot_z");
  if (parit != m_GlobalParameterMap.end())
  {
    _rot_in_z = parit->second * cm;
  }

  return 0;
}
