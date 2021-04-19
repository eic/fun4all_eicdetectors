#include "PHG4ForwardDualReadoutDetector.h"
#include "PHG4ForwardDualReadoutDisplayAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
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
  , _place_in_x(0.0 * mm)
  , _place_in_y(0.0 * mm)
  , _place_in_z(4000.0 * mm)
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
  , _tower_dx(100 * mm)
  , _tower_dy(100 * mm)
  , _tower_dz(1000.0 * mm)
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
    /* only record energy in actual absorber- drop energy lost in air gaps inside drcalo envelope */
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

  /* Read parameters for detector construction and mappign from file */
  InitDefaultParams();

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* drcalo_envelope_solid = new G4Cons("hdrcalo_envelope_solid",
                                             _rMin1, _rMax1,
                                             _rMin2, _rMax2,
                                             _dZ / 2.0,
                                             _sPhi, _dPhi);

  G4LogicalVolume* drcalo_envelope_log = new G4LogicalVolume(drcalo_envelope_solid, Air, G4String("hdrcalo_envelope"), 0, 0, 0);

  m_DisplayAction->AddVolume(drcalo_envelope_log, "FdrcaloEnvelope");

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix drcalo_rotm;
  drcalo_rotm.rotateX(_rot_in_x);
  drcalo_rotm.rotateY(_rot_in_y);
  drcalo_rotm.rotateZ(_rot_in_z);

  /* Place envelope cone in simulation */
  ostringstream name_envelope;
  name_envelope.str("");
  name_envelope << _towerlogicnameprefix << "_envelope" << endl;

  new G4PVPlacement(G4Transform3D(drcalo_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z)),
                    drcalo_envelope_log, name_envelope.str().c_str(), logicWorld, 0, false, OverlapCheck());


  G4LogicalVolume* singletower = ConstructTower(0); //4x4 fibre tower with 2 scint and 2 Cherenkov fibres


  G4Material* material_air = G4Material::GetMaterial("G4_AIR");
 // number of towers in radial direction (on y axis)
  int rowNtow = (int) ( (_rMax1-(_tower_dy/2)) / _tower_dy);
  for(int row=rowNtow;row>=-rowNtow;row--){
    // pythagoras -> get available length in circular mother volume for towers
    // divide given length by tower width -> get number of towers that can be placed
    int currRowNtow = (int) ( ( 2* sqrt(pow(_rMax1,2)-pow( (abs(row)*_tower_dy) + (_tower_dy/2.0) ,2)) ) / _tower_dy );
    // we want an odd number of towers to be symmetrically centered around 0
    if ( currRowNtow % 2 != 0) currRowNtow-=1;

    // TODO account for hole at center of detector
    // float rinner = 50.0 * mm;
    if( ( (abs(row)*_tower_dy) - (_tower_dy/2.0) ) < _rMin1 ){ // _rMin1

      // pythagoras -> get available length in circular mother volume for towers
      // divide given length by tower width -> get number of towers that can be placed
      int currRowNtowInner = (int) ( ( 2* sqrt(pow(_rMin1,2)-pow( (abs(row)*_tower_dy) - (_tower_dy/2.0) ,2)) ) / _tower_dy );
      // we want an odd number of towers to be symmetrically centered around 0
      if ( currRowNtowInner % 2 != 0) currRowNtowInner-=1;
      // currRowNtowInner+=2;

      // create mother volume with space for currRowNtow towers along x-axis
      auto HadCalRowSolid    = new G4Box("HadCalRowBox", (currRowNtow - currRowNtowInner) / 2 * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
      auto HadCalRowLogical  = new G4LogicalVolume(HadCalRowSolid,material_air,"HadCalRowLogical");
      // replicate singletower tower design currRowNtow times along x-axis
      new G4PVReplica("HadCalRowPhysical",singletower,HadCalRowLogical,
                      kXAxis,(currRowNtow - currRowNtowInner) / 2,_tower_dx);

      ostringstream name_row_twr;
      name_row_twr.str("");
      name_row_twr << _towerlogicnameprefix << "_row_" << row << "_left" << endl;
      new G4PVPlacement(0, G4ThreeVector( - ( ( currRowNtowInner / 2.0 ) * _tower_dx ) - ( (currRowNtow - currRowNtowInner) / 2 * _tower_dx / 2.0 ), (row*_tower_dy), 0),
                    HadCalRowLogical, name_row_twr.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());

      ostringstream name_row_twr2;
      name_row_twr2.str("");
      name_row_twr2 << _towerlogicnameprefix << "_row_" << row << "_left" << endl;
      new G4PVPlacement(0, G4ThreeVector( ( ( currRowNtowInner / 2.0 ) * _tower_dx ) + ( (currRowNtow - currRowNtowInner) / 2 * _tower_dx / 2.0 ), (row*_tower_dy), 0),
                    HadCalRowLogical, name_row_twr2.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());

    } else {
      // create mother volume with space for currRowNtow towers along x-axis
      auto HadCalRowSolid    = new G4Box("HadCalRowBox", currRowNtow * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
      auto HadCalRowLogical  = new G4LogicalVolume(HadCalRowSolid,material_air,"HadCalRowLogical");
      // replicate singletower tower design currRowNtow times along x-axis
      new G4PVReplica("HadCalRowPhysical",singletower,HadCalRowLogical,
                      kXAxis,currRowNtow,_tower_dx);
      
      ostringstream name_row_twr;
      name_row_twr.str("");
      name_row_twr << _towerlogicnameprefix << "_row_" << row << endl;
      new G4PVPlacement(0, G4ThreeVector(0, (row*_tower_dy), 0),
                    HadCalRowLogical, name_row_twr.str().c_str(), drcalo_envelope_log, 0, false, OverlapCheck());
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

  // if (type == 1) return ConstructTowerType2();

  /* create logical volume for single tower */
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");
  float distancing = 1.0;
  // 2x2 tower base element
  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"),
                                          _tower_dx / 2.0,
                                          _tower_dy / 2.0,
                                          _tower_dz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  /* create geometry volumes to place inside single_tower */

  G4double diameter_fiber = 1 * mm;
  G4double airgap = 0.1 * mm;

  // fibre cutout
  G4VSolid* single_cutout_tube  = new G4Tubs(G4String("single_cutout_tube"),
                                            0,
                                            ( diameter_fiber + airgap ) / 2.0,
                                            1.03 * _tower_dz / 1.0, //make it 1.03 times longer to ensure full cutout
                                            0.,2*M_PI*rad);
  // notch cutout
  G4VSolid* single_cutout_box = new G4Box(G4String("single_cutout_box"),
                                          ( diameter_fiber + airgap ) / 2.0,
                                          1.03 * ( diameter_fiber + airgap ) / 4.0, //make it 1.03 times longer to ensure full cutout
                                          1.03 * _tower_dz / 1.0);
  // absorber base object
  G4VSolid* solid_absorber_temp = new G4Box(G4String("solid_absorber_temp"),
                                          distancing*_tower_dx / 2.0,
                                          distancing*_tower_dy / 2.0,
                                          _tower_dz / 2.0);
  // cut out four fiber holes
  G4VSolid* solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_temp_f1"), solid_absorber_temp, single_cutout_tube
                                                            , 0 ,G4ThreeVector( _tower_dx / 4.0 , ( _tower_dy / 2.0 ) - ( ( diameter_fiber + airgap ) / 2.0 ) ,0.)); // top right
  solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_temp_f2"), solid_absorber, single_cutout_tube
                                                            , 0 ,G4ThreeVector( _tower_dx / 4.0 , - ( ( diameter_fiber + airgap ) / 2.0 ) ,0.)); // bottom right
  solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_temp_f3"), solid_absorber, single_cutout_tube
                                                            , 0 ,G4ThreeVector(- _tower_dx / 4.0, ( _tower_dy / 2.0 ) - ( ( diameter_fiber + airgap ) / 2.0 ) ,0.)); // top left
  solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_temp_f4"), solid_absorber, single_cutout_tube
                                                            , 0 ,G4ThreeVector(- _tower_dx / 4.0, - ( ( diameter_fiber + airgap ) / 2.0 ) ,0.)); // bottom left
  // G4VSolid* solid_absorber = solid_absorber_fiber->Clone();
  // cut out four notches
  solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_temp_box1"), solid_absorber, single_cutout_box
                                                            , 0 ,G4ThreeVector( _tower_dx / 4.0  , ( _tower_dy / 2.0 ) - ( ( diameter_fiber + airgap ) / 4.0 ) ,0.));
  solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_temp_box2"), solid_absorber, single_cutout_box
                                                            , 0 ,G4ThreeVector( _tower_dx / 4.0  , - ( ( diameter_fiber + airgap ) / 4.0 ) ,0.));
  solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_temp_box3"), solid_absorber, single_cutout_box
                                                            , 0 ,G4ThreeVector( - _tower_dx / 4.0, ( _tower_dy / 2.0 ) - ( ( diameter_fiber + airgap ) / 4.0 ) ,0.));
  solid_absorber = new G4SubtractionSolid(G4String("solid_absorber_temp_box4"), solid_absorber, single_cutout_box
                                                            , 0 ,G4ThreeVector( - _tower_dx / 4.0, - ( ( diameter_fiber + airgap ) / 4.0 ) ,0.));



  G4VSolid* solid_scintillator  = new G4Tubs(G4String("single_scintillator_fiber"),
                                            0,
                                            diameter_fiber / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PI*rad);


  G4VSolid* solid_cherenkov  = new G4Tubs(G4String("single_cherenkov_fiber"),
                                            0,
                                            diameter_fiber / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PI*rad);



// const G4int NUMENTRIES = 9;
// G4double LXe_PP[NUMENTRIES] = {6.6*eV,6.7*eV,6.8*eV,6.9*eV,7.0*eV,7.1*eV,7.2*eV,7.3*eV,7.4*eV};
// G4double LXe_SCINT[NUMENTRIES] = {0.000134, 0.004432, 0.053991,0.241971, 0.398942, 0.000134, 0.004432, 0.053991,0.241971};
// G4double LXe_RIND[NUMENTRIES] = { 1.57, 1.57, 1.57, 1.57, 1.57, 1.57,1.57, 1.57, 1.57};
// G4double LXe_ABSL[NUMENTRIES] = { 35.*cm, 35.*cm, 35.*cm, 35.*cm,35.*cm, 35.*cm, 35.*cm, 35.*cm, 35*cm };

// G4MaterialPropertiesTable* LXe_MPT = new G4MaterialPropertiesTable();
// LXe_MPT -> AddProperty("FASTCOMPONENT",LXe_PP,LXe_SCINT,NUMENTRIES);
// LXe_MPT -> AddProperty("RINDEX",LXe_PP,LXe_RIND,NUMENTRIES);
// LXe_MPT -> AddProperty("ABSLENGTH",LXe_PP,LXe_ABSL,NUMENTRIES);
// LXe_MPT -> AddConstProperty ("SCINTILLATIONYIELD",10000./MeV);
// LXe_MPT -> AddConstProperty("RESOLUTIONSCALE",1.0);
// LXe_MPT -> AddConstProperty("FASTTIMECONSTANT",45.*ns);
// LXe_MPT -> AddConstProperty("YIELDRATIO",1.0);

  // const G4int ntab = 2;
  // G4double scin_en[] = {2.9*eV, 3.*eV}; // 420 nm (the range is 414 - 428 nm)
  // G4double scin_fast[] = {1., 1.};

  G4MaterialPropertiesTable *tab = new G4MaterialPropertiesTable();

  if(1){
    const G4int ntab = 31;
    // tab->AddProperty("FASTCOMPONENT", scin_en, scin_fast, ntab);
    tab->AddConstProperty("FASTTIMECONSTANT", 2.8*ns); // was 6
    // tab->AddConstProperty("SCINTILLATIONYIELD", 13.9/keV); // was 200/MEV nominal  10
    tab->AddConstProperty("SCINTILLATIONYIELD", 200/MeV); // was 200/MEV nominal  10
    tab->AddConstProperty("RESOLUTIONSCALE", 1.);

    // G4double opt_en[] = {1.551*eV, 3.545*eV}; // 350 - 800 nm
    // G4double opt_r[] = {2.4, 2.4};
    // G4double opt_abs[] = {200*cm, 200*cm};

    // tab->AddProperty("RINDEX", opt_en, opt_r, ntab);
    // tab->AddProperty("ABSLENGTH", opt_en, opt_abs, ntab);
    G4double opt_en[] = {1.37760*eV, 1.45864*eV, 1.54980*eV, 1.65312*eV, 1.71013*eV, 1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV, 2.10143*eV, 2.13766*eV, 2.17516*eV, 2.21400*eV, 2.25426*eV, 2.29600*eV, 2.33932*eV, 2.38431*eV, 2.43106*eV, 2.47968*eV, 2.53029*eV, 2.58300*eV, 2.63796*eV, 2.69531*eV, 2.75520*eV, 2.81782*eV, 2.88335*eV, 2.95200*eV, 3.09960*eV, 3.54241*eV, 4.13281*eV}; // 350 - 800 nm
    G4double scin_fast[] = {0, 0, 0, 0, 0, 0, 0, 0.0003, 0.0008, 0.0032, 0.0057, 0.0084, 0.0153, 0.0234, 0.0343, 0.0604, 0.0927, 0.1398, 0.2105, 0.2903, 0.4122, 0.5518, 0.7086, 0.8678, 1, 0.8676, 0.2311, 0.0033, 0.0012, 0, 0};
    tab->AddProperty("FASTCOMPONENT", opt_en, scin_fast, ntab);

    G4double opt_r[] = {1.5749, 1.5764, 1.5782, 1.5803, 1.5815, 1.5829, 1.5845, 1.5862, 1.5882, 1.5904, 1.5914, 1.5924, 1.5935, 1.5947, 1.5959, 1.5972, 1.5986, 1.6, 1.6016, 1.6033, 1.6051, 1.607, 1.609, 1.6112, 1.6136, 1.6161, 1.617, 1.623, 1.62858, 1.65191, 1.69165};
    tab->AddProperty("RINDEX", opt_en, opt_r, ntab);

    G4double opt_abs[] = {2.714*m, 3.619*m, 5.791*m, 4.343*m, 7.896*m, 5.429*m, 36.19*m, 17.37*m, 36.19*m, 5.429*m, 13.00*m, 14.50*m, 16.00*m, 18.00*m, 16.50*m, 17.00*m, 14.00*m, 16.00*m, 15.00*m, 14.50*m, 13.00*m, 12.00*m, 10.00*m, 8.000*m, 7.238*m, 4.000*m, 1.200*m, 0.500*m, 0.200*m, 0.200*m, 0.100*m};
    tab->AddProperty("ABSLENGTH", opt_en, opt_abs, ntab);

      //   <material name="DR_Polystyrene">
      // <D value="1.032" unit="g/cm3"/>
      // <composite n="19" ref="C"/>
  }

  G4Material* material_scintillator = G4Material::GetMaterial(_materialScintillator.c_str());
  material_scintillator->SetMaterialPropertiesTable(tab);
  // material_scintillator->SetMaterialPropertiesTable(LXe_MPT);

  // G4Material* material_absorber = G4Material::GetMaterial(_materialAbsorber.c_str());
  G4NistManager* man = G4NistManager::Instance();
  G4Material* material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());


  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "absorber_solid_logic",
                                                        0, 0, 0);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                    material_scintillator,
                                                    "hdrcalo_single_scintillator_fiber_logic",
                                                    0, 0, 0);



  // G4OpticalSurface *surface = new G4OpticalSurface("CrystalSurface", unified, polished, dielectric_metal);
  // if(1){
  //   //G4LogicalSkinSurface *csurf = 
  //   new G4LogicalSkinSurface("CrystalSurfaceL", logic_scint, surface);

  //   //surface material
  //   const G4int ntab = 2;
  //   G4double opt_en[] = {1.551*eV, 3.545*eV}; // 350 - 800 nm
  //   G4double reflectivity[] = {0.8, 0.8};
  //   G4double efficiency[] = {0.9, 0.9};
  //   G4MaterialPropertiesTable *surfmat = new G4MaterialPropertiesTable();
  //   surfmat->AddProperty("REFLECTIVITY", opt_en, reflectivity, ntab);
  //   surfmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  //   surface->SetMaterialPropertiesTable(surfmat);
  // }

  // G4Material* material_cherenkov = G4Material::GetMaterial("PMMA");
	G4double density;
	G4int ncomponents, natoms;
  G4Material* material_cherenkov = new G4Material("PMMA", density = 1.18 * g / cm3, ncomponents = 3);
		material_cherenkov->AddElement(G4Element::GetElement("C"), 3.6 / (3.6 + 5.7 + 1.4));
		material_cherenkov->AddElement(G4Element::GetElement("H"), 5.7 / (3.6 + 5.7 + 1.4));
		material_cherenkov->AddElement(G4Element::GetElement("O"), 1.4 / (3.6 + 5.7 + 1.4));


const G4int nEntries = 31;

  G4double photonEnergy[nEntries] =
      {1.37760*eV, 1.45864*eV, 1.54980*eV, 1.65312*eV, 1.71013*eV, 1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV, 2.10143*eV, 2.13766*eV, 2.17516*eV, 2.21400*eV, 2.25426*eV, 2.29600*eV, 2.33932*eV, 2.38431*eV, 2.43106*eV, 2.47968*eV, 2.53029*eV, 2.58300*eV, 2.63796*eV, 2.69531*eV, 2.75520*eV, 2.81782*eV, 2.88335*eV, 2.95200*eV, 3.09960*eV, 3.54241*eV, 4.13281*eV};
    G4double refractiveIndexWLSfiber[nEntries] =
    { 1.4852, 1.4859, 1.4867, 1.4877, 1.4882, 1.4888, 1.4895, 1.4903, 1.4911, 1.492, 1.4924, 1.4929, 1.4933, 1.4938, 1.4943, 1.4948, 1.4954, 1.496, 1.4966, 1.4973, 1.4981, 1.4989, 1.4997, 1.5006, 1.5016, 1.5026, 1.5038, 1.505, 1.5052, 1.5152, 1.5306};

  G4double absWLSfiber[nEntries] =
      {0.414*m, 0.965*m, 2.171*m, 4.343*m, 1.448*m, 4.343*m, 14.48*m, 21.71*m, 8.686*m, 39.48*m, 48.25*m, 54.29*m, 57.91*m, 54.29*m, 33.40*m, 31.02*m, 43.43*m, 43.43*m, 41.36*m, 39.48*m, 37.76*m, 36.19*m, 36.19*m, 33.40*m, 31.02*m, 28.95*m, 25.55*m, 24.13*m, 21.71*m, 2.171*m, 0.434*m};


  // Add entries into properties table
  G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
  mptWLSfiber->AddProperty("RINDEX",photonEnergy,refractiveIndexWLSfiber,nEntries);
  mptWLSfiber->AddProperty("ABSLENGTH",photonEnergy,absWLSfiber,nEntries);
  material_cherenkov->SetMaterialPropertiesTable(mptWLSfiber);


  G4LogicalVolume* logic_cherenk = new G4LogicalVolume(solid_cherenkov,
                                                    material_cherenkov,
                                                    "hdrcalo_single_cherenkov_fiber_logic",
                                                    0, 0, 0);

  m_DisplayAction->AddVolume(logic_absorber, "Absorber");
  m_DisplayAction->AddVolume(logic_scint, "Scintillator");
  m_DisplayAction->AddVolume(logic_cherenk, "Cherenkov");

  /* place physical volumes for absorber and scintillator fiber */

  // ostringstream towername;
  // towername.str("");
  // towername ;

  ostringstream name_absorber;
  name_absorber.str("");
  name_absorber << _towerlogicnameprefix << "absorbersolid" << endl;

  ostringstream name_scintillator;
  name_scintillator.str("");
  name_scintillator << _towerlogicnameprefix << "singlescintillatorfiber" << endl;

  ostringstream name_cherenkov;
  name_cherenkov.str("");
  name_cherenkov << _towerlogicnameprefix << "_single_cherenkov_fiber"  << endl;


  new G4PVPlacement(0, G4ThreeVector( 0,  0 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  // place scintillator fibers (top left, bottom right)
  new G4PVPlacement(0, G4ThreeVector( -_tower_dx / 4,  ( _tower_dy / 2.0 ) - ( ( diameter_fiber + airgap ) / 2.0 ) , 0),
                    logic_scint,
                    name_scintillator.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( _tower_dx / 4,  - ( ( diameter_fiber + airgap ) / 2.0 ) , 0),
                    logic_scint,
                    name_scintillator.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  // place cherenkov fibers (top right, bottom left)
  new G4PVPlacement(0, G4ThreeVector( _tower_dx / 4 ,  ( _tower_dy / 2.0 ) - ( ( diameter_fiber + airgap ) / 2.0 ) , 0),
                    logic_cherenk,
                    name_cherenkov.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( -_tower_dx / 4 ,  - ( ( diameter_fiber + airgap ) / 2.0 ) , 0),
                    logic_cherenk,
                    name_cherenkov.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());



  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardDualReadoutDetector: Building logical volume for single tower done." << endl;
  }

  return single_tower_logic;
}


int PHG4ForwardDualReadoutDetector::InitDefaultParams()
{
    _tower_dx = 0.3 * cm;
    _tower_dy = 0.3 * cm;
    _tower_dz = 150 * cm; //was 150
    _rMin1 = 20 * cm;
    _rMax1 = 220 * cm; //was220
    _rMin2 = 20 * cm;
    _rMax2 = 220 * cm; //was220
    _dZ = 150 * cm; //was 150
    _place_in_x = 0 * cm;
    _place_in_y = 0 * cm;
    _place_in_z = 375 * cm;
    _rot_in_x = 0;
    _rot_in_y = 0;
    _rot_in_z = 0;

  return 0;
}