#include "PHG4BSTDetector.h"
#include "PHG4BSTDisplayAction.h"
#include "PHG4BSTSteppingAction.h"

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
PHG4BSTDetector::PHG4BSTDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4BSTDisplayAction*>(subsys->GetDisplayAction()))
  , m_SteppingAction(0)
  , _place_in_x(0.0 * mm)
  , _place_in_y(0.0 * mm)
  , _place_in_z(0.0 * mm)
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
  , _towerlogicnameprefix("hbstTower")
  , _superdetector("NONE")
  , _mapping_tower_file("")
{
}
//_______________________________________________________________________
int PHG4BSTDetector::IsInForwardDualReadout(G4VPhysicalVolume* volume) const
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
    //only record energy in actual absorber- drop energy lost in air gaps inside bst envelope
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
void PHG4BSTDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4BSTDetector: Begin Construction" << endl;
  }

  //Read parameters for detector construction from file
  // ParseParametersFromTable();

  //Create the cone envelope = 'world volume' for the calorimeter
  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4double minR = 3.0 * cm;
  G4double maxR = 29.0* cm;
  G4double lengthBST = 90.0* cm;
  G4VSolid* bst_envelope_solid = new G4Cons("hbst_envelope_solid",
                                        minR, maxR,
                                        minR, maxR,
                                        lengthBST / 2.0,
                                        _sPhi, _dPhi);

  G4LogicalVolume* bst_envelope_log = new G4LogicalVolume(bst_envelope_solid, Air, G4String("hbst_envelope"), 0, 0, 0);

  m_DisplayAction->AddVolume(bst_envelope_log, "Invisible");
  // m_DisplayAction->AddVolume(bst_envelope_log, "FbstEnvelope");

  //Define rotation attributes for envelope cone
  G4RotationMatrix bst_rotm;
  bst_rotm.rotateX(_rot_in_x);
  bst_rotm.rotateY(_rot_in_y);
  bst_rotm.rotateZ(_rot_in_z);

  //Place envelope cone in simulation
  ostringstream name_envelope;
  name_envelope.str("");
  name_envelope << _towerlogicnameprefix << "_envelope" << endl;

  new G4PVPlacement(G4Transform3D(bst_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z)),
                    bst_envelope_log, name_envelope.str().c_str(), logicWorld, 0, false, OverlapCheck());

  ConstructBarrel(bst_envelope_log);
  // ConstructOuterBarrel(bst_envelope_log);
  return;
}

void PHG4BSTDetector::ConstructBarrel(G4LogicalVolume* mother){


  G4double layer_sensor_thickness = 0.05 / 100 * 9.37 * cm; // 0.05% of Si rad length 9.37cm
  G4Material *layer_material = GetDetectorMaterial("G4_Si");
  G4double deadarea_seam = 0.1 * cm;

  const int nLayersInner = 3;
  G4double layer_radius_inner[nLayersInner] = {3.40 * cm, 5.67 * cm, 7.93 * cm};
  G4double layer_length_inner[nLayersInner] = {30.0 * cm, 30.0 * cm, 30.0 * cm};


  const int nLayersOuter = 2;
  G4double layer_radius_outer[nLayersOuter] = {15.3 * cm, 17.0 * cm};
  G4double layer_length_outer[nLayersOuter] = {30.0 * cm, 30.0 * cm}; // two layers of 30cm
  int layer_segments[nLayersOuter] = {12, 12};
  G4double copperWire_diam = 0.64 * mm;

  G4double support_radius_inner = 9.0 * cm;
  G4double support_length_inner = 2 * 19.44 * cm;
  G4double support_thickness_foam = 0.2 * cm;
  G4double support_thickness_shell = 0.1 * mm;
  G4double support_thickness_shell_outer = 0.2 * mm;
  // G4double support_seam = 0.1 * cm;

  G4double support_radius_outer = 18.5 * cm;
  G4double support_length_outer = 2 * 37.7861 * cm;

  G4Material *foam_material = MakeCarbonFoamMaterial();
  G4double foam_length = 1.0 * cm;
  G4double foam_spacing = 1.3 * cm;
  G4double foam_thickness = 0.25 * cm;
  G4double foam_endwheel_depth = 1.3 * cm;
  int foam_endwheel_holes_inner = 12;
  G4double foam_endwheel_hole_diam_inner[nLayersInner] = {(layer_radius_inner[1]-layer_radius_inner[0])/3, (layer_radius_inner[2]-layer_radius_inner[1])/3, (support_radius_inner-layer_radius_inner[2])/2};
  int foam_endwheel_holes_outer = 30;
  G4double foam_endwheel_hole_diam_outer[nLayersOuter] = {(layer_radius_outer[1]-layer_radius_outer[0])/3, (support_radius_outer-layer_radius_outer[1])/3};


  for(int i = 0; i < nLayersInner; i++){
    G4double deadangle_seam = deadarea_seam / layer_radius_inner[i];
    G4VSolid* currentLayerSolid  = new G4Tubs(G4String("currentVertexLayerSolid"),
                                            layer_radius_inner[i] - layer_sensor_thickness / 2,
                                            layer_radius_inner[i] + layer_sensor_thickness / 2,
                                            layer_length_inner[i] / 2,
                                            0.+deadangle_seam,M_PI*rad-deadangle_seam);
  
    G4LogicalVolume* currentLayerLogic = new G4LogicalVolume(currentLayerSolid,
                                                          layer_material,
                                                          "currentVertexLayerLogic_"+std::to_string(i),
                                                          0, 0, 0);
 
    m_DisplayAction->AddVolume(currentLayerLogic, "InnerBarrel");


    G4double angle_foam = foam_thickness / layer_radius_inner[i];

    G4VSolid* foamblockSolid  = new G4Tubs(G4String("foamblockSolid"),
                                            layer_radius_inner[i] + layer_sensor_thickness / 2,
                                            i<nLayersInner-1 ? layer_radius_inner[i+1] - layer_sensor_thickness / 2 : support_radius_inner - support_thickness_foam / 2 - support_thickness_shell,
                                            foam_length / 2,
                                            -angle_foam/2,angle_foam/2);
  
    G4LogicalVolume* foamblockLogic = new G4LogicalVolume(foamblockSolid,
                                                          foam_material,
                                                          "foamblockLogic_"+std::to_string(i),
                                                          0, 0, 0);
 
    m_DisplayAction->AddVolume(foamblockLogic, "Foam");

    G4VSolid* foamEndWheelSolid  = new G4Tubs(G4String("foamEndWheelSolid"),
                                            layer_radius_inner[i] + layer_sensor_thickness / 2,
                                            i<nLayersInner-1 ? layer_radius_inner[i+1] - layer_sensor_thickness / 2 : support_radius_inner - support_thickness_foam / 2 - support_thickness_shell,
                                            foam_endwheel_depth / 2,
                                            0.,M_PI*rad);
                                            // 0.+deadangle_seam,M_PI*rad-deadangle_seam);
    G4VSolid* foamEndWheelStencilSolid  = new G4Tubs(G4String("foamEndWheelStencilSolid"),
                                            0,
                                            foam_endwheel_hole_diam_inner[i] / 2,
                                            foam_endwheel_depth,
                                            0.,2*M_PI);
    for(int j = 1; j < foam_endwheel_holes_inner; j++){
      foamEndWheelSolid = new G4SubtractionSolid(G4String("foamEndWheelSolid"+std::to_string(j)),
                                            foamEndWheelSolid,
                                            foamEndWheelStencilSolid, 0,
                                            G4ThreeVector((layer_radius_inner[i]+(i<nLayersInner-1 ? layer_radius_inner[i+1] : support_radius_inner))/2 * cos(M_PI*j/foam_endwheel_holes_inner),(layer_radius_inner[i]+(i<nLayersInner-1 ? layer_radius_inner[i+1] : support_radius_inner))/2 * sin(M_PI*j/foam_endwheel_holes_inner),0));
    }
    G4LogicalVolume* foamEndWheelLogic = new G4LogicalVolume(foamEndWheelSolid,
                                                          foam_material,
                                                          "foamEndWheelLogic_"+std::to_string(i),
                                                          0, 0, 0);
 
    m_DisplayAction->AddVolume(foamEndWheelLogic, "FoamEndWheel");

    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, 0.0),
                      currentLayerLogic,
                      "currentVertexLayerLogicTop"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());
    G4RotationMatrix bstlayer_rotm;
    bstlayer_rotm.rotateZ(M_PI);

    new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, 0.0)),
                      currentLayerLogic,
                      "currentVertexLayerLogicBottom"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());

    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, -layer_length_inner[i] / 2 + foam_length / 2),
                      foamEndWheelLogic,
                      "foamEndWheelPlacedFrontTop"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, layer_length_inner[i] / 2 - foam_length / 2),
                      foamEndWheelLogic,
                      "foamEndWheelPlacedBackTop"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());


    new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, -layer_length_inner[i] / 2 + foam_length / 2)),
                      foamEndWheelLogic,
                      "foamEndWheelPlacedFrontBottom"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, layer_length_inner[i] / 2 - foam_length / 2)),
                      foamEndWheelLogic,
                      "foamEndWheelPlacedBackBottom"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());

    double addRotate[3] = {1.5*angle_foam, 0, -1.5*angle_foam};
    for(int ifoamphi=0; ifoamphi<3; ifoamphi++){
      for(int ifoamz=1; ifoamz<(int)(layer_length_inner[i]/(foam_length+foam_spacing))-1; ifoamz++){
        G4RotationMatrix foam_rotm;
        foam_rotm.rotateZ(ifoamphi*M_PI/2 + addRotate[ifoamphi]+deadangle_seam/2);
        new G4PVPlacement(G4Transform3D(foam_rotm, G4ThreeVector( 0.0, 0.0, -layer_length_inner[i]/2 + (foam_length+foam_spacing)*(ifoamz+0.5))),
                          foamblockLogic,
                          "foamblockLogicTop_"+std::to_string(ifoamz)+"_"+std::to_string(ifoamphi),
                          mother,
                          0, 0, OverlapCheck());

        G4RotationMatrix foam_rotm2;
        foam_rotm2.rotateZ(ifoamphi*M_PI/2 + addRotate[ifoamphi]+M_PI+deadangle_seam/2);
        new G4PVPlacement(G4Transform3D(foam_rotm2, G4ThreeVector( 0.0, 0.0, -layer_length_inner[i]/2 + (foam_length+foam_spacing)*(ifoamz+0.5))),
                          foamblockLogic,
                          "foamblockLogicBottom_"+std::to_string(ifoamz)+"_"+std::to_string(ifoamphi),
                          mother,
                          0, 0, OverlapCheck());
      }

    }
  }

  for(int i = 0; i < nLayersOuter; i++){
    G4double deadangle_seam = deadarea_seam / layer_radius_outer[i];
    G4VSolid* currentLayerSolid  = new G4Tubs("currentSagittaLayerSolid"+std::to_string(i),
                                            layer_radius_outer[i] - layer_sensor_thickness / 2,
                                            layer_radius_outer[i] + layer_sensor_thickness / 2,
                                            layer_length_outer[i] / 2,
                                            0.+deadangle_seam,(2*M_PI*rad/layer_segments[i])-deadangle_seam);
  
    G4LogicalVolume* currentLayerLogic = new G4LogicalVolume(currentLayerSolid,
                                                          layer_material,
                                                          "currentSagittaLayerLogic_"+std::to_string(i),
                                                          0, 0, 0);

    m_DisplayAction->AddVolume(currentLayerLogic, "OuterBarrel");

    G4VSolid* copperWireSolid  = new G4Tubs("copperWireSolid"+std::to_string(i),
                                            0,
                                            copperWire_diam / 2,
                                            layer_length_outer[i] / 2 -  support_thickness_foam,
                                            0.,(2*M_PI*rad));
  
    G4LogicalVolume* copperWireLogic = new G4LogicalVolume(copperWireSolid,
                                                          GetDetectorMaterial("G4_Cu", false),
                                                          "copperWireLogic_"+std::to_string(i),
                                                          0, 0, 0);

    m_DisplayAction->AddVolume(copperWireLogic, "CopperWire");
    for(int iseg = 0; iseg < layer_segments[i]; iseg++){
      G4RotationMatrix bstlayer_rotm;
      bstlayer_rotm.rotateZ(2*M_PI*iseg/layer_segments[i]);

      new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, deadarea_seam+ layer_length_outer[i] / 2)),
                        currentLayerLogic,
                        "currentSagittaLayerLogicFront"+std::to_string(iseg),
                        mother,
                        0, 0, OverlapCheck());
      new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, -(deadarea_seam+layer_length_outer[i] / 2))),
                        currentLayerLogic,
                        "currentSagittaLayerLogicBack"+std::to_string(iseg),
                        mother,
                        0, 0, OverlapCheck());
      if((iseg-1)%2==0){
        new G4PVPlacement(0, G4ThreeVector((layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * cos(2*M_PI*iseg/layer_segments[i]),(layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * sin(2*M_PI*iseg/layer_segments[i]),layer_length_outer[i] / 2),
                          copperWireLogic,
                          "copperWireLogicBack"+std::to_string(iseg),
                          mother,
                          0, 0, OverlapCheck());
        new G4PVPlacement(0, G4ThreeVector((layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * cos(2*M_PI*iseg/layer_segments[i]),(layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * sin(2*M_PI*iseg/layer_segments[i]),-layer_length_outer[i] / 2),
                          copperWireLogic,
                          "copperWireLogicFront"+std::to_string(iseg),
                          mother,
                          0, 0, OverlapCheck());
      }
    }

    G4VSolid* foamEndWheelSagittaSolid  = new G4Tubs(G4String("foamEndWheelSagittaSolid"),
                                            layer_radius_outer[i] + layer_sensor_thickness / 2,
                                            i<nLayersOuter-1 ? layer_radius_outer[i+1] - layer_sensor_thickness / 2 : support_radius_outer - support_thickness_foam / 2 - support_thickness_shell_outer,
                                            foam_endwheel_depth / 2,
                                            0.,M_PI*rad);
                                            // 0.+deadangle_seam,M_PI*rad-deadangle_seam);
    G4VSolid* foamEndWheelSagittaStencilSolid  = new G4Tubs(G4String("foamEndWheelSagittaStencilSolid"),
                                            0,
                                            foam_endwheel_hole_diam_outer[i] / 2,
                                            foam_endwheel_depth,
                                            0.,2*M_PI);
    for(int j = 1; j < foam_endwheel_holes_outer; j++){
      foamEndWheelSagittaSolid = new G4SubtractionSolid(G4String("foamEndWheelSagittaSolid"+std::to_string(j)),
                                            foamEndWheelSagittaSolid,
                                            foamEndWheelSagittaStencilSolid, 0,
                                            G4ThreeVector((layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * cos(M_PI*j/foam_endwheel_holes_outer),(layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * sin(M_PI*j/foam_endwheel_holes_outer),0));
    }
    G4LogicalVolume* foamEndWheelSagittaLogic = new G4LogicalVolume(foamEndWheelSagittaSolid,
                                                          foam_material,
                                                          "foamEndWheelSagittaLogic_"+std::to_string(i),
                                                          0, 0, 0);
    m_DisplayAction->AddVolume(foamEndWheelSagittaLogic, "FoamEndWheel");

    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, -layer_length_outer[i] + foam_length / 2),
                      foamEndWheelSagittaLogic,
                      "foamEndWheelSagittaPlacedFrontTop"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, 0.0),
                      foamEndWheelSagittaLogic,
                      "foamEndWheelSagittaPlacedCenterTop"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, layer_length_outer[i] - foam_length / 2),
                      foamEndWheelSagittaLogic,
                      "foamEndWheelSagittaPlacedBackTop"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());

    G4RotationMatrix bstlayer_rotm180;
    bstlayer_rotm180.rotateZ(M_PI);
    new G4PVPlacement(G4Transform3D(bstlayer_rotm180, G4ThreeVector( 0.0, 0.0, -layer_length_outer[i] + foam_length / 2)),
                      foamEndWheelSagittaLogic,
                      "foamEndWheelSagittaPlacedFrontBottom"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(G4Transform3D(bstlayer_rotm180, G4ThreeVector( 0.0, 0.0, 0.0)),
                      foamEndWheelSagittaLogic,
                      "foamEndWheelSagittaPlacedCenterBottom"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(G4Transform3D(bstlayer_rotm180, G4ThreeVector( 0.0, 0.0, layer_length_outer[i] - foam_length / 2)),
                      foamEndWheelSagittaLogic,
                      "foamEndWheelSagittaPlacedBackBottom"+std::to_string(i),
                      mother,
                      0, 0, OverlapCheck());

  }

  G4Material *support_material = GetCarbonFiber();
  G4double support_seamangle_inner = 0;//support_seam / support_radius_inner;

  // construct support 2mm thick foam with 0.1mm carbon skins
  G4VSolid* supportCylinderFoamSolid  = new G4Tubs(G4String("supportVertexCylinderSolid"),
                                          support_radius_inner - support_thickness_foam / 2,
                                          support_radius_inner + support_thickness_foam / 2,
                                          support_length_inner / 2,
                                          0.+support_seamangle_inner,M_PI*rad-support_seamangle_inner);
  G4LogicalVolume* supportCylinderFoamLogic = new G4LogicalVolume(supportCylinderFoamSolid,
                                                        foam_material,
                                                        "supportCylinderFoamLogic",
                                                        0, 0, 0);

  G4VSolid* supportCylinderOuterShellSolid  = new G4Tubs(G4String("supportVertexCylinderSolid"),
                                          support_radius_inner + support_thickness_foam / 2,
                                          support_radius_inner + support_thickness_foam / 2 + support_thickness_shell,
                                          support_length_inner / 2,
                                          0.+support_seamangle_inner,M_PI*rad-support_seamangle_inner);
  G4LogicalVolume* supportCylinderOuterShellLogic = new G4LogicalVolume(supportCylinderOuterShellSolid,
                                                        support_material,
                                                        "supportCylinderOuterShellLogic",
                                                        0, 0, 0);
 
  G4VSolid* supportCylinderInnerShellSolid  = new G4Tubs(G4String("supportVertexCylinderSolid"),
                                          support_radius_inner - support_thickness_foam / 2 - support_thickness_shell,
                                          support_radius_inner - support_thickness_foam / 2,
                                          support_length_inner / 2,
                                          0.+support_seamangle_inner,M_PI*rad-support_seamangle_inner);
  G4LogicalVolume* supportCylinderInnerShellLogic = new G4LogicalVolume(supportCylinderInnerShellSolid,
                                                        support_material,
                                                        "supportCylinderInnerShellLogic",
                                                        0, 0, 0);
 
    m_DisplayAction->AddVolume(supportCylinderFoamLogic, "Foam");
    m_DisplayAction->AddVolume(supportCylinderOuterShellLogic, "CShell");
    m_DisplayAction->AddVolume(supportCylinderInnerShellLogic, "CShell");


    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, 0.0),
                      supportCylinderFoamLogic,
                      "supportCylinderFoamLogicLogicTop",
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, 0.0),
                      supportCylinderOuterShellLogic,
                      "supportCylinderOuterShellLogicLogicTop",
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, 0.0),
                      supportCylinderInnerShellLogic,
                      "supportCylinderInnerShellLogicLogicTop",
                      mother,
                      0, 0, OverlapCheck());
    G4RotationMatrix bstsupp_rotm;
    bstsupp_rotm.rotateZ(M_PI);
    new G4PVPlacement(G4Transform3D(bstsupp_rotm, G4ThreeVector( 0.0, 0.0, 0.0)),
                      supportCylinderFoamLogic,
                      "supportCylinderFoamLogicLogicBottom",
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(G4Transform3D(bstsupp_rotm, G4ThreeVector( 0.0, 0.0, 0.0)),
                      supportCylinderOuterShellLogic,
                      "supportCylinderOuterShellLogicLogicBottom",
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(G4Transform3D(bstsupp_rotm, G4ThreeVector( 0.0, 0.0, 0.0)),
                      supportCylinderInnerShellLogic,
                      "supportCylinderInnerShellLogicLogicBottom",
                      mother,
                      0, 0, OverlapCheck());


  G4double support_seamangle_outer = 0;//support_seam / support_radius_outer;

  // construct support 2mm thick foam with 0.1mm carbon skins
  G4VSolid* supportSagittaFoamSolid  = new G4Tubs(G4String("supportVertexCylinderSolid"),
                                          support_radius_outer - support_thickness_foam / 2,
                                          support_radius_outer + support_thickness_foam / 2,
                                          support_length_outer / 2,
                                          0.+support_seamangle_outer,M_PI*rad-support_seamangle_outer);
  G4LogicalVolume* supportSagittaFoamLogic = new G4LogicalVolume(supportSagittaFoamSolid,
                                                        foam_material,
                                                        "supportSagittaFoamLogic",
                                                        0, 0, 0);

  G4VSolid* supportSagittaOuterShellSolid  = new G4Tubs(G4String("supportVertexCylinderSolid"),
                                          support_radius_outer + support_thickness_foam / 2,
                                          support_radius_outer + support_thickness_foam / 2 + support_thickness_shell_outer,
                                          support_length_outer / 2,
                                          0.+support_seamangle_outer,M_PI*rad-support_seamangle_outer);
  G4LogicalVolume* supportSagittaOuterShellLogic = new G4LogicalVolume(supportSagittaOuterShellSolid,
                                                        support_material,
                                                        "supportSagittaOuterShellLogic",
                                                        0, 0, 0);
 
  G4VSolid* supportSagittaInnerShellSolid  = new G4Tubs(G4String("supportVertexCylinderSolid"),
                                          support_radius_outer - support_thickness_foam / 2 - support_thickness_shell_outer,
                                          support_radius_outer - support_thickness_foam / 2,
                                          support_length_outer / 2,
                                          0.+support_seamangle_outer,M_PI*rad-support_seamangle_outer);
  G4LogicalVolume* supportSagittaInnerShellLogic = new G4LogicalVolume(supportSagittaInnerShellSolid,
                                                        support_material,
                                                        "supportSagittaInnerShellLogic",
                                                        0, 0, 0);
 
    m_DisplayAction->AddVolume(supportSagittaFoamLogic, "Foam");
    m_DisplayAction->AddVolume(supportSagittaOuterShellLogic, "CShell");
    m_DisplayAction->AddVolume(supportSagittaInnerShellLogic, "CShell");


    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, 0.0),
                      supportSagittaFoamLogic,
                      "supportSagittaFoamLogicLogicTop",
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, 0.0),
                      supportSagittaOuterShellLogic,
                      "supportSagittaOuterShellLogicLogicTop",
                      mother,
                      0, 0, OverlapCheck());
                      
    new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, 0.0),
                      supportSagittaInnerShellLogic,
                      "supportSagittaInnerShellLogicLogicTop",
                      mother,
                      0, 0, OverlapCheck());


    new G4PVPlacement(G4Transform3D(bstsupp_rotm, G4ThreeVector( 0.0, 0.0, 0.0)),
                      supportSagittaFoamLogic,
                      "supportSagittaFoamLogicLogicBottom",
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(G4Transform3D(bstsupp_rotm, G4ThreeVector( 0.0, 0.0, 0.0)),
                      supportSagittaOuterShellLogic,
                      "supportSagittaOuterShellLogicLogicBottom",
                      mother,
                      0, 0, OverlapCheck());
    new G4PVPlacement(G4Transform3D(bstsupp_rotm, G4ThreeVector( 0.0, 0.0, 0.0)),
                      supportSagittaInnerShellLogic,
                      "supportSagittaInnerShellLogicLogicBottom",
                      mother,
                      0, 0, OverlapCheck());

  return;
}

G4Material* PHG4BSTDetector::MakeCarbonFoamMaterial(){
  G4Material* carbon_foam = GetDetectorMaterial("C_FOAM_BST", false);  // false suppresses warning that material does not exist
  if(!carbon_foam){
    G4double density;
    G4int ncomponents;
    carbon_foam = new G4Material("C_FOAM_BST", density = 0.26 * g / cm3, ncomponents = 2);
    // carbon_foam = new G4Material("C_FOAM_BST", density = 0.06 * g / cm3, ncomponents = 2);
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.97);
    carbon_foam->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 0.03);
  }
  return carbon_foam;

}

G4Material* PHG4BSTDetector::GetCarbonFiber()
{
  static string matname = "BSTCarbonFiber";
  G4Material* carbonfiber = G4Material::GetMaterial(matname, false);  // false suppresses warning that material does not exist
  if (!carbonfiber)
  {
    G4double density_carbon_fiber = 1.44 * g / cm3;
    carbonfiber = new G4Material(matname, density_carbon_fiber, 1);
    carbonfiber->AddElement(G4Element::GetElement("C"), 1);
  }
  return carbonfiber;
}

int PHG4BSTDetector::ParseParametersFromTable()
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
