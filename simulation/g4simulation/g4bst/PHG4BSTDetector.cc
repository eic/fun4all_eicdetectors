#include "PHG4BSTDetector.h"
#include "PHG4BSTDisplayAction.h"
#include "PHG4BSTSteppingAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phparameter/PHParameters.h>

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
PHG4BSTDetector::PHG4BSTDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters *parameters, const std::string& dnam)
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
  , m_Params(parameters)
{
}
//_______________________________________________________________________
int PHG4BSTDetector::IsInActiveSensorBST(G4VPhysicalVolume* volume) const
{
  // if (volume->GetName().find(_towerlogicnameprefix) != string::npos)
  // {
    if (volume->GetName().find("currentVertexLayer") != string::npos)
    {
      if (_active)
        return 1;
      else
        return 0;
    }
    else if (volume->GetName().find("currentSagittaLayer") != string::npos)
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
  // }

  return 0;
}

//_______________________________________________________________________
void PHG4BSTDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4BSTDetector: Begin Construction" << endl;
  }

  ConstructBarrel(logicWorld);
  return;
}

void PHG4BSTDetector::ConstructBarrel(G4LogicalVolume* mother){

  bool do_internal_supports = true;
  bool do_external_supports = true;

  // Sensor sizes available from the ALPIDE wafer
  G4double sensor_widths[4] = {56.5 * mm, 75.5 * mm, 94.3 * mm, 113.1 *mm};
  G4double sensor_length[4] = {270 * mm, 270 * mm, 270 * mm, 240 *mm};
  G4double layer_sensor_thickness = 0.05 / 100 * 9.37 * cm; // 0.05% of Si rad length 9.37cm // 46.85 microns
  G4Material *layer_material = GetDetectorMaterial("G4_Si");
  G4double layer_backing_thickness = m_Params->get_double_param("layer_backing_thickness"); //0.07 / 100 * 28.57 * cm; // 200 microns of Kapton
  G4Material *backing_material = GetDetectorMaterial("G4_KAPTON");

  // default seam where two half shells connect
  G4double deadarea_seam = 0.5 * mm;
  // default seam where two half shells connect
  G4double deadarea_sensoredge = 0.5 * mm;


  G4double support_radius_inner = 7.5 * cm; // was 9 cm
  G4double support_length_inner = 38.88 * cm;
  G4double support_thickness_foam = 0.2 * cm;
  G4double support_thickness_shell = 0.1 * mm;
  G4double support_thickness_shell_outer = 0.2 * mm;
  // G4double support_seam = 0.1 * cm;

  G4double support_radius_outer = 22.5 * cm; // was 18.5
  G4double support_length_outer = 2 * 37.7861 * cm; // was 2 * 37.7861

  // 3 vertex layers, 2 sagitta layers
  const int nLayersInner = 3;
  G4double layer_radius_inner[nLayersInner] = {0};
  G4double layer_length_inner[nLayersInner] = {0};
  int nSensorsInner[nLayersInner] = {0};
  G4double layer_sensor_width_inner[nLayersInner][9] = {0};

  const int nLayersOuter = 2;
  G4double layer_radius_outer[nLayersOuter] = {0};
  G4double layer_length_outer[nLayersOuter] = {0};
  int nSensorsOuter[nLayersOuter] = {0};
  G4double layer_sensor_width_outer[nLayersOuter][9] = {0};

  // vertex layers sensor setup
  nSensorsInner[0] = 1;
  layer_radius_inner[0] = (sensor_widths[3]+deadarea_seam)/M_PI; // 36.319158 mm -> eta 1.92
  layer_sensor_width_inner[0][0] = sensor_widths[3];
  layer_length_inner[0] = sensor_length[3];

  bool use_different_vertex_sensors = false;
  if(use_different_vertex_sensors){
    nSensorsInner[0] = 2;
    layer_radius_inner[0] = (2*sensor_widths[0]+deadarea_seam)/M_PI; // 36.287327 mm -> eta 2.0
    layer_sensor_width_inner[0][0] = sensor_widths[0];
    layer_sensor_width_inner[0][1] = sensor_widths[0];
    layer_length_inner[0] = sensor_length[0];
  }

  nSensorsInner[1] = 2;
  layer_radius_inner[1] = (2*sensor_widths[1]+deadarea_seam)/M_PI; // 48.383103 mm
  layer_sensor_width_inner[1][0] = sensor_widths[1];
  layer_sensor_width_inner[1][1] = sensor_widths[1];
  layer_length_inner[1] = sensor_length[1];

  nSensorsInner[2] = 2;
  layer_radius_inner[2] = (2*sensor_widths[2]+deadarea_seam)/M_PI; // 60.351554 mm -> eta 1.49
  layer_sensor_width_inner[2][0] = sensor_widths[2];
  layer_sensor_width_inner[2][1] = sensor_widths[2];
  layer_length_inner[2] = sensor_length[2];

  // sagitta layers sensor setup
  // Sagitta1(nonproj+): 5*94.3+2*75.5 = radius 198 -> eta 0.95
  // Sagitta2(nonproj+): 7*94.3 = radius 210 -> eta 0.89
  nSensorsOuter[0] = 7;
  layer_radius_outer[0] = (5*sensor_widths[2]+2*sensor_widths[1]+deadarea_seam)/M_PI; // 198.46621 mm
  for(int i=0;i<5;i++){
    layer_sensor_width_outer[0][i] = sensor_widths[2];
  }
  layer_sensor_width_outer[0][5] = sensor_widths[1];
  layer_sensor_width_outer[0][6] = sensor_widths[1];
  layer_length_outer[0] = sensor_length[2];

  nSensorsOuter[1] = 7;
  layer_radius_outer[1] = (7*sensor_widths[2]+deadarea_seam)/M_PI; // 210.43467 mm
  for(int i=0;i<7;i++){
    layer_sensor_width_outer[1][i] = sensor_widths[2];
  }
  layer_length_outer[1] = sensor_length[2];


  bool use_sagitta_setup_1 = false;
  if(use_sagitta_setup_1){
    // Sagitta1: 5*94.3 = radius 150.09 -> eta 1.25
    // Sagitta2: 94.3*4+75.5*2 = radius 168.13 -> eta 1.13
    // OR 94.3*5+56.5 = radius 168.07 -> eta 1.13
    nSensorsOuter[0] = 5;
    layer_radius_outer[0] = (5*sensor_widths[2]+deadarea_seam)/M_PI; // 150.40142 mm
    for(int i=0;i<5;i++){
      layer_sensor_width_outer[0][i] = sensor_widths[2];
    }
    layer_length_outer[0] = sensor_length[2];


    nSensorsOuter[1] = 6;
    layer_radius_outer[1] = (4*sensor_widths[2]+2*sensor_widths[1]+deadarea_seam)/M_PI; // 168.44959 mm
    for(int i=0;i<4;i++){
      layer_sensor_width_outer[1][i] = sensor_widths[2];
    }
    layer_sensor_width_outer[1][4] = sensor_widths[1];
    layer_sensor_width_outer[1][5] = sensor_widths[1];
    layer_length_outer[1] = sensor_length[2];

    support_radius_outer = 18.2 * cm;
  }

  bool use_sagitta_setup_2 = false;
  if(use_sagitta_setup_2){
    // Sagitta1(nonproj): 7*94.3 = radius 210 -> eta 0.89
    // Sagitta2(nonproj): 6*94.3+2*75.5 = radius 228 -> eta 0.77
    nSensorsOuter[0] = 7;
    layer_radius_outer[0] = (7*sensor_widths[2]+deadarea_seam)/M_PI; // 210.43467 mm
    for(int i=0;i<7;i++){
      layer_sensor_width_outer[0][i] = sensor_widths[2];
    }
    layer_length_outer[0] = sensor_length[2];


    nSensorsOuter[1] = 8;
    layer_radius_outer[1] = (6*sensor_widths[2]+2*sensor_widths[1]+deadarea_seam)/M_PI; // 228.48284 mm
    for(int i=0;i<6;i++){
      layer_sensor_width_outer[1][i] = sensor_widths[2];
    }
    layer_sensor_width_outer[1][6] = sensor_widths[1];
    layer_sensor_width_outer[1][7] = sensor_widths[1];
    layer_length_outer[1] = sensor_length[2];

    support_radius_outer = 24.2 * cm;
  }
  for(int i=0; i<nLayersInner; i++){
    cout << "BST layer " << i << " radius " << layer_radius_inner[i] << " length " << layer_length_inner[i] << " nSensors " << nSensorsInner[i] << endl;
    // for(int j=0; j<nSensorsInner[i]; j++){
    //   cout << "\tlayer " << i << " sensor " << j << " width " << layer_sensor_width_inner[i][j] << endl;
    // }
  }
  for(int i=0; i<nLayersOuter; i++){
    cout << "BST layer " << i+nLayersInner << " radius " << layer_radius_outer[i] << " length " << layer_length_outer[i] << " nSensors " << nSensorsOuter[i] << endl;
    // for(int j=0; j<nSensorsOuter[i]; j++){
    //   cout << "\tlayer " << i+nLayersInner << " sensor " << j << " width " << layer_sensor_width_outer[i][j] << endl;
    // }
  }
  // int layer_segments[nLayersOuter] = {12, 12};
  G4double copperWire_diam = 0.64 * mm;


  G4Material *foam_material = MakeCarbonFoamMaterial();
  G4double foam_length = 0.8 * cm;
  // G4double foam_spacing = 1.3 * cm;
  G4double foam_thickness = 0.4 * cm;
  G4double foam_endwheel_depth = 1.0 * cm;
  G4double foam_midwheel_depth = 0.7 * cm;
  int foam_endwheel_holes_inner = 16;
  G4double foam_endwheel_hole_diam_inner[nLayersInner] = {(layer_radius_inner[1]-layer_radius_inner[0])/2.5, (layer_radius_inner[2]-layer_radius_inner[1])/2.5, (support_radius_inner-layer_radius_inner[2])/2};
  int foam_endwheel_holes_outer = 60;
  G4double foam_endwheel_hole_diam_outer[nLayersOuter] = {(layer_radius_outer[1]-layer_radius_outer[0])/2, (support_radius_outer-layer_radius_outer[1])/2.5};


  for(int i = 0; i < nLayersInner; i++){
    G4double layer_half_circumference = M_PI * layer_radius_inner[i];
    G4double angle_deadarea_sensoredge = M_PI * (deadarea_sensoredge / 2.0) / layer_half_circumference;

    G4double current_angle = M_PI * (deadarea_seam / 2.0) / layer_half_circumference;
    current_angle += angle_deadarea_sensoredge;

    for(int j = 0; j < nSensorsInner[i]; j++){
      G4double angle_current_layer = M_PI * (layer_sensor_width_inner[i][j] - deadarea_sensoredge) / layer_half_circumference;
      // cout << "layer_sensor_width_inner[" << i << "][" << j << "] = " << layer_sensor_width_inner[i][j] << "\tangle_current_layer = " << angle_current_layer << "\tstarting_angle = " << current_angle << endl;
      G4VSolid* currentLayerSolid  = new G4Tubs(G4String( "BST_"+std::to_string(i)+"_currentVertexLayerSolid_" + std::to_string(i) + "_" + std::to_string(j)),
                                            layer_radius_inner[i] - layer_sensor_thickness / 2,
                                            layer_radius_inner[i] + layer_sensor_thickness / 2,
                                            layer_length_inner[i] / 2,
                                            0,angle_current_layer);

      G4LogicalVolume* currentLayerLogic = new G4LogicalVolume(currentLayerSolid,
                                                            layer_material,
                                                            "BST_"+std::to_string(i)+"_currentVertexLayerLogic_"+std::to_string(i) + "_" + std::to_string(j),
                                                            0, 0, 0);
  
      m_DisplayAction->AddVolume(currentLayerLogic, "InnerBarrel");

      G4RotationMatrix bstlayer_rotm;
      bstlayer_rotm.rotateZ(current_angle);
      new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, 0.0)),
                        currentLayerLogic,
                        "BST_"+std::to_string(i)+"_currentVertexLayerPlacedTop"+std::to_string(i)+ "_" + std::to_string(j),
                        mother,
                        0, 0, OverlapCheck());

      G4RotationMatrix bstlayer_rotm_bottom;
      bstlayer_rotm_bottom.rotateZ(current_angle+M_PI);

      new G4PVPlacement(G4Transform3D(bstlayer_rotm_bottom, G4ThreeVector( 0.0, 0.0, 0.0)),
                        currentLayerLogic,
                        "BST_"+std::to_string(i)+"_currentVertexLayerPlacedBottom"+std::to_string(i)+ "_" + std::to_string(j),
                        mother,
                        0, 0, OverlapCheck());
      current_angle+=(angle_current_layer+2*angle_deadarea_sensoredge);


      if(do_internal_supports){
        G4double angle_foam = M_PI * foam_thickness / layer_half_circumference;

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


        double addRotate[3] = {-1.5*angle_foam, 0, 1.5*angle_foam};
        for(int ifoamphi=0; ifoamphi<3; ifoamphi++){
          G4RotationMatrix foam_rotm;
          if(ifoamphi==0){
            foam_rotm.rotateZ(current_angle - (angle_deadarea_sensoredge) + addRotate[ifoamphi]);
          } else if(ifoamphi==1){
            foam_rotm.rotateZ(current_angle - (angle_deadarea_sensoredge + angle_current_layer/2.0) + addRotate[ifoamphi]);
          } else if(ifoamphi==2){
            foam_rotm.rotateZ(current_angle - (angle_deadarea_sensoredge + angle_current_layer) + addRotate[ifoamphi]);
          }
          G4RotationMatrix foam_rotm2;
          if(ifoamphi==0){
            foam_rotm2.rotateZ(M_PI+current_angle - (angle_deadarea_sensoredge) + addRotate[ifoamphi]);
          } else if(ifoamphi==1){
            foam_rotm2.rotateZ(M_PI+current_angle - (angle_deadarea_sensoredge + angle_current_layer/2.0) + addRotate[ifoamphi]);
          } else if(ifoamphi==2){
            foam_rotm2.rotateZ(M_PI+current_angle - (angle_deadarea_sensoredge + angle_current_layer) + addRotate[ifoamphi]);
          }
          G4double spaceForSpacers = layer_length_inner[0] / 3.0;
          for(int ispacerow=-1; ispacerow<2; ispacerow++){
            for(int ifoamz=-1; ifoamz<2; ifoamz++){
              new G4PVPlacement(G4Transform3D(foam_rotm, G4ThreeVector( 0.0, 0.0, ispacerow*spaceForSpacers + ifoamz*(spaceForSpacers/3.0))),
                                foamblockLogic,
                                "foamblockLogicTop_"+std::to_string(ispacerow)+"_"+std::to_string(ifoamz)+"_"+std::to_string(ifoamphi)+"_"+std::to_string(i)+ "_" + std::to_string(j),
                                mother,
                                0, 0, OverlapCheck());
              new G4PVPlacement(G4Transform3D(foam_rotm2, G4ThreeVector( 0.0, 0.0, ispacerow*spaceForSpacers + ifoamz*(spaceForSpacers/3.0))),
                                foamblockLogic,
                                "foamblockLogicBottom_"+std::to_string(ispacerow)+"_"+std::to_string(ifoamz)+"_"+std::to_string(ifoamphi)+"_"+std::to_string(i)+ "_" + std::to_string(j),
                                mother,
                                0, 0, OverlapCheck());

              // new G4PVPlacement(G4Transform3D(foam_rotm2, G4ThreeVector( 0.0, 0.0, -layer_length_inner[i]/2 + (foam_length+foam_spacing)*(ifoamz+0.5))),
              //                   foamblockLogic,
              //                   "foamblockLogicBottom_"+std::to_string(ifoamz)+"_"+std::to_string(ifoamphi)+"_"+std::to_string(i)+ "_" + std::to_string(j),
              //                   mother,
              //                   0, 0, OverlapCheck());
            }
          }
        }
      }
    }


    if(do_internal_supports){
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

      G4VSolid* foamMidWheelSolid  = new G4Tubs(G4String("foamMidWheelSolid"),
                                              layer_radius_inner[i] + layer_sensor_thickness / 2,
                                              i<nLayersInner-1 ? layer_radius_inner[i+1] - layer_sensor_thickness / 2 : support_radius_inner - support_thickness_foam / 2 - support_thickness_shell,
                                              foam_midwheel_depth / 2,
                                              0.,M_PI*rad);
                                              // 0.+deadangle_seam,M_PI*rad-deadangle_seam);
      G4VSolid* foamMidWheelStencilSolid  = new G4Tubs(G4String("foamMidWheelStencilSolid"),
                                              0,
                                              foam_endwheel_hole_diam_inner[i] / 2,
                                              foam_midwheel_depth,
                                              0.,2*M_PI);
      for(int j = 1; j < foam_endwheel_holes_inner; j++){
        foamMidWheelSolid = new G4SubtractionSolid(G4String("foamMidWheelSolid"+std::to_string(j)),
                                              foamMidWheelSolid,
                                              foamMidWheelStencilSolid, 0,
                                              G4ThreeVector((layer_radius_inner[i]+(i<nLayersInner-1 ? layer_radius_inner[i+1] : support_radius_inner))/2 * cos(M_PI*j/foam_endwheel_holes_inner),(layer_radius_inner[i]+(i<nLayersInner-1 ? layer_radius_inner[i+1] : support_radius_inner))/2 * sin(M_PI*j/foam_endwheel_holes_inner),0));
      }
      G4LogicalVolume* foamMidWheelLogic = new G4LogicalVolume(foamMidWheelSolid,
                                                            foam_material,
                                                            "foamMidWheelLogic_"+std::to_string(i),
                                                            0, 0, 0);
  
      m_DisplayAction->AddVolume(foamMidWheelLogic, "FoamEndWheel");


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

      new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, -layer_length_inner[0] / 6 + foam_length / 2),
                        foamMidWheelLogic,
                        "foamMidWheelPlacedFrontTop"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, layer_length_inner[0] / 6 - foam_length / 2),
                        foamMidWheelLogic,
                        "foamMidWheelPlacedBackTop"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());

      G4RotationMatrix bstlayer_rotm_wheel;
      bstlayer_rotm_wheel.rotateZ(M_PI);

      new G4PVPlacement(G4Transform3D(bstlayer_rotm_wheel, G4ThreeVector( 0.0, 0.0, -layer_length_inner[i] / 2 + foam_length / 2)),
                        foamEndWheelLogic,
                        "foamEndWheelPlacedFrontBottom"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());
      new G4PVPlacement(G4Transform3D(bstlayer_rotm_wheel, G4ThreeVector( 0.0, 0.0, layer_length_inner[i] / 2 - foam_length / 2)),
                        foamEndWheelLogic,
                        "foamEndWheelPlacedBackBottom"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());

      new G4PVPlacement(G4Transform3D(bstlayer_rotm_wheel, G4ThreeVector( 0.0, 0.0, -layer_length_inner[0] / 6 + foam_length / 2)),
                        foamMidWheelLogic,
                        "foamMidWheelPlacedFrontBottom"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());
      new G4PVPlacement(G4Transform3D(bstlayer_rotm_wheel, G4ThreeVector( 0.0, 0.0, layer_length_inner[0] / 6 - foam_length / 2)),
                        foamMidWheelLogic,
                        "foamMidWheelPlacedBackBottom"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());
    }
  }

  for(int i = 0; i < nLayersOuter; i++){
    G4double layer_half_circumference = M_PI * layer_radius_outer[i];
    G4double angle_deadarea_sensoredge = M_PI * (deadarea_sensoredge / 2.0) / layer_half_circumference;

    G4double current_angle = M_PI * (deadarea_seam / 2.0) / layer_half_circumference;
    current_angle += angle_deadarea_sensoredge;

    for(int j = 0; j < nSensorsOuter[i]; j++){
      G4double angle_current_layer = M_PI * (layer_sensor_width_outer[i][j] - deadarea_sensoredge) / layer_half_circumference;
      // cout << "layer_sensor_width_outer[" << i << "][" << j << "] = " << layer_sensor_width_outer[i][j] << "\tangle_current_layer = " << angle_current_layer << "\tstarting_angle = " << current_angle << endl;
      G4VSolid* currentLayerSolid  = new G4Tubs(G4String("BST_"+std::to_string(nLayersInner+i)+"_currentSagittaLayerSolid_" + std::to_string(i) + "_" + std::to_string(j)),
                                            layer_radius_outer[i] - layer_sensor_thickness / 2,
                                            layer_radius_outer[i] + layer_sensor_thickness / 2,
                                            layer_length_outer[i] / 2,
                                            0,angle_current_layer);

      G4LogicalVolume* currentLayerLogic = new G4LogicalVolume(currentLayerSolid,
                                                            layer_material,
                                                            "BST_"+std::to_string(nLayersInner+i)+"_currentSagittaLayerLogic_"+std::to_string(i) + "_" + std::to_string(j),
                                                            0, 0, 0);
  
      m_DisplayAction->AddVolume(currentLayerLogic, "InnerBarrel");

      G4LogicalVolume* currentLayerBackingLogic = nullptr;
      if(layer_backing_thickness>0){
        G4VSolid* currentLayerBackingSolid  = new G4Tubs(G4String("SagittaBackingLayerSolid_" + std::to_string(i) + "_" + std::to_string(j)),
                                              layer_radius_outer[i] + layer_sensor_thickness / 2,
                                              layer_radius_outer[i] + layer_sensor_thickness / 2 + layer_backing_thickness,
                                              layer_length_outer[i] / 2,
                                              0,angle_current_layer);

        currentLayerBackingLogic = new G4LogicalVolume(currentLayerBackingSolid,
                                                              backing_material,
                                                              "SagittaBackingLayerLogic_"+std::to_string(i) + "_" + std::to_string(j),
                                                              0, 0, 0);
    
        m_DisplayAction->AddVolume(currentLayerBackingLogic, "BackingMaterial");
      }


      G4RotationMatrix bstlayer_rotm;
      bstlayer_rotm.rotateZ(current_angle);
      new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, deadarea_seam / 2.0 + layer_length_outer[i] / 2)),
                        currentLayerLogic,
                        "BST_"+std::to_string(nLayersInner+i)+"_currentSagittaLayerPlacedTopBack_"+std::to_string(i)+ "_" + std::to_string(j),
                        mother,
                        0, 0, OverlapCheck());
      new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, -(deadarea_seam / 2.0 + layer_length_outer[i] / 2))),
                        currentLayerLogic,
                        "BST_"+std::to_string(nLayersInner+i)+"_currentSagittaLayerPlacedTopFront_"+std::to_string(i)+ "_" + std::to_string(j),
                        mother,
                        0, 0, OverlapCheck());
      if(currentLayerBackingLogic){
        new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, deadarea_seam / 2.0 + layer_length_outer[i] / 2)),
                          currentLayerBackingLogic,
                          "BackingSagittaLayerPlacedTopBack_"+std::to_string(i)+ "_" + std::to_string(j),
                          mother,
                          0, 0, OverlapCheck());
        new G4PVPlacement(G4Transform3D(bstlayer_rotm, G4ThreeVector( 0.0, 0.0, -(deadarea_seam / 2.0 + layer_length_outer[i] / 2))),
                          currentLayerBackingLogic,
                          "BackingSagittaLayerPlacedTopFront_"+std::to_string(i)+ "_" + std::to_string(j),
                          mother,
                          0, 0, OverlapCheck());
      }
      G4RotationMatrix bstlayer_rotm_bottom;
      bstlayer_rotm_bottom.rotateZ(current_angle+M_PI);

      new G4PVPlacement(G4Transform3D(bstlayer_rotm_bottom, G4ThreeVector( 0.0, 0.0, deadarea_seam / 2.0 + layer_length_outer[i] / 2)),
                        currentLayerLogic,
                        "BST_"+std::to_string(nLayersInner+i)+"_currentSagittaLayerPlacedBottomBack_"+std::to_string(i)+ "_" + std::to_string(j),
                        mother,
                        0, 0, OverlapCheck());
      new G4PVPlacement(G4Transform3D(bstlayer_rotm_bottom, G4ThreeVector( 0.0, 0.0, -(deadarea_seam / 2.0 + layer_length_outer[i] / 2))),
                        currentLayerLogic,
                        "BST_"+std::to_string(nLayersInner+i)+"_currentSagittaLayerPlacedBottomFront_"+std::to_string(i)+ "_" + std::to_string(j),
                        mother,
                        0, 0, OverlapCheck());
      if(currentLayerBackingLogic){
        new G4PVPlacement(G4Transform3D(bstlayer_rotm_bottom, G4ThreeVector( 0.0, 0.0, deadarea_seam / 2.0 + layer_length_outer[i] / 2)),
                          currentLayerBackingLogic,
                          "BackingSagittaLayerPlacedBottomBack_"+std::to_string(i)+ "_" + std::to_string(j),
                          mother,
                          0, 0, OverlapCheck());
        new G4PVPlacement(G4Transform3D(bstlayer_rotm_bottom, G4ThreeVector( 0.0, 0.0, -(deadarea_seam / 2.0 + layer_length_outer[i] / 2))),
                          currentLayerBackingLogic,
                          "BackingSagittaLayerPlacedBottomFront_"+std::to_string(i)+ "_" + std::to_string(j),
                          mother,
                          0, 0, OverlapCheck());
      }

      G4VSolid* copperWireSolid  = new G4Tubs("copperWireSolid_"+std::to_string(i)+ "_" + std::to_string(j),
                                              0,
                                              copperWire_diam / 2,
                                              layer_length_outer[i] / 4 -  foam_endwheel_depth,
                                              0.,(2*M_PI*rad));
    
      G4LogicalVolume* copperWireLogic = new G4LogicalVolume(copperWireSolid,
                                                            GetDetectorMaterial("G4_Cu", false),
                                                            "copperWireLogic_"+std::to_string(i)+ "_" + std::to_string(j),
                                                            0, 0, 0);

      m_DisplayAction->AddVolume(copperWireLogic, "CopperWire");

      if((j)%2==0){
        // place copper wire between every second sensor
        for(int k=-2; k<3; k++){
          if(k==0) continue;
          new G4PVPlacement(0, G4ThreeVector((layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * cos(current_angle),(layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * sin(current_angle),k*layer_length_outer[i] / 2 + (k<0 ? layer_length_outer[i]/4 : -layer_length_outer[i]/4)),
                                    copperWireLogic,
                                    "copperWireLogicTop_"+std::to_string(i)+ "_" + std::to_string(j)+ "_" + std::to_string(k),
                                    mother,
                                    0, 0, OverlapCheck());
          new G4PVPlacement(0, G4ThreeVector((layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * cos(current_angle+M_PI),(layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * sin(current_angle+M_PI),k*layer_length_outer[i] / 2 + (k<0 ? layer_length_outer[i]/4 : -layer_length_outer[i]/4)),
                                    copperWireLogic,
                                    "copperWireLogicBottom_"+std::to_string(i)+ "_" + std::to_string(j)+ "_" + std::to_string(k),
                                    mother,
                                    0, 0, OverlapCheck());
        }
      }


      if(do_internal_supports){
        G4double angle_foam = M_PI * foam_thickness / layer_half_circumference;

        G4VSolid* foamblockSolid  = new G4Tubs(G4String("foamblockSolid"),
                                                layer_radius_outer[i] + layer_sensor_thickness / 2 + layer_backing_thickness,
                                                i<nLayersOuter-1 ? layer_radius_outer[i+1] - layer_sensor_thickness / 2 : support_radius_outer - support_thickness_foam / 2 - support_thickness_shell_outer,
                                                foam_length / 2,
                                                -angle_foam/2,angle_foam/2);
      
        G4LogicalVolume* foamblockLogic = new G4LogicalVolume(foamblockSolid,
                                                              foam_material,
                                                              "foamblockSagittaLogic_"+std::to_string(i),
                                                              0, 0, 0);
    
        m_DisplayAction->AddVolume(foamblockLogic, "Foam");

        current_angle+=(angle_current_layer+2*angle_deadarea_sensoredge);

        double addRotate[3] = {-1.5*angle_foam, 0, 1.5*angle_foam};
        for(int ifoamphi=0; ifoamphi<3; ifoamphi++){
          G4RotationMatrix foam_rotm;
          if(ifoamphi==0){
            foam_rotm.rotateZ(current_angle - (angle_deadarea_sensoredge) + addRotate[ifoamphi]);
          } else if(ifoamphi==1){
            foam_rotm.rotateZ(current_angle - (angle_deadarea_sensoredge + angle_current_layer/2.0) + addRotate[ifoamphi]);
          } else if(ifoamphi==2){
            foam_rotm.rotateZ(current_angle - (angle_deadarea_sensoredge + angle_current_layer) + addRotate[ifoamphi]);
          }
          G4RotationMatrix foam_rotm2;
          if(ifoamphi==0){
            foam_rotm2.rotateZ(M_PI+current_angle - (angle_deadarea_sensoredge) + addRotate[ifoamphi]);
          } else if(ifoamphi==1){
            foam_rotm2.rotateZ(M_PI+current_angle - (angle_deadarea_sensoredge + angle_current_layer/2.0) + addRotate[ifoamphi]);
          } else if(ifoamphi==2){
            foam_rotm2.rotateZ(M_PI+current_angle - (angle_deadarea_sensoredge + angle_current_layer) + addRotate[ifoamphi]);
          }
          G4double spaceForSpacers = layer_length_outer[i] / 2.0;
          for(int ispacerow=0; ispacerow<4; ispacerow++){
            for(int ifoamz=-2; ifoamz<3; ifoamz++){
              new G4PVPlacement(G4Transform3D(foam_rotm, G4ThreeVector( 0.0, 0.0, - layer_length_outer[i] + spaceForSpacers/2.0 + ispacerow*spaceForSpacers + ifoamz*(spaceForSpacers/5.0))),
                                foamblockLogic,
                                "foamblockSagittaPlacedTop_"+std::to_string(ispacerow)+"_"+std::to_string(ifoamz)+"_"+std::to_string(ifoamphi)+"_"+std::to_string(i)+ "_" + std::to_string(j),
                                mother,
                                0, 0, OverlapCheck());
              new G4PVPlacement(G4Transform3D(foam_rotm2, G4ThreeVector( 0.0, 0.0, - layer_length_outer[i] + spaceForSpacers/2.0 + ispacerow*spaceForSpacers + ifoamz*(spaceForSpacers/5.0))),
                                foamblockLogic,
                                "foamblockSagittaPlacedBottom_"+std::to_string(ispacerow)+"_"+std::to_string(ifoamz)+"_"+std::to_string(ifoamphi)+"_"+std::to_string(i)+ "_" + std::to_string(j),
                                mother,
                                0, 0, OverlapCheck());

              // new G4PVPlacement(G4Transform3D(foam_rotm2, G4ThreeVector( 0.0, 0.0, -layer_length_inner[i]/2 + (foam_length+foam_spacing)*(ifoamz+0.5))),
              //                   foamblockLogic,
              //                   "foamblockLogicBottom_"+std::to_string(ifoamz)+"_"+std::to_string(ifoamphi)+"_"+std::to_string(i)+ "_" + std::to_string(j),
              //                   mother,
              //                   0, 0, OverlapCheck());
            }
          }
        }
      }
    }

    if(do_internal_supports){
      G4VSolid* foamEndWheelSagittaSolid  = new G4Tubs(G4String("foamEndWheelSagittaSolid"),
                                              layer_radius_outer[i] + layer_sensor_thickness / 2 + layer_backing_thickness,
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

      G4VSolid* foamMidWheelSagittaSolid  = new G4Tubs(G4String("foamMidWheelSagittaSolid"),
                                              layer_radius_outer[i] + layer_sensor_thickness / 2 + layer_backing_thickness,
                                              i<nLayersOuter-1 ? layer_radius_outer[i+1] - layer_sensor_thickness / 2 : support_radius_outer - support_thickness_foam / 2 - support_thickness_shell_outer,
                                              foam_midwheel_depth / 2,
                                              0.,M_PI*rad);
                                              // 0.+deadangle_seam,M_PI*rad-deadangle_seam);
      G4VSolid* foamMidWheelSagittaStencilSolid  = new G4Tubs(G4String("foamMidWheelSagittaStencilSolid"),
                                              0,
                                              foam_endwheel_hole_diam_outer[i] / 2,
                                              foam_endwheel_depth,
                                              0.,2*M_PI);
      for(int j = 1; j < foam_endwheel_holes_outer; j++){
        foamMidWheelSagittaSolid = new G4SubtractionSolid(G4String("foamMidWheelSagittaSolid"+std::to_string(j)),
                                              foamMidWheelSagittaSolid,
                                              foamMidWheelSagittaStencilSolid, 0,
                                              G4ThreeVector((layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * cos(M_PI*j/foam_endwheel_holes_outer),(layer_radius_outer[i]+(i<nLayersOuter-1 ? layer_radius_outer[i+1] : support_radius_outer))/2 * sin(M_PI*j/foam_endwheel_holes_outer),0));
      }
      G4LogicalVolume* foamMidWheelSagittaLogic = new G4LogicalVolume(foamMidWheelSagittaSolid,
                                                            foam_material,
                                                            "foamMidWheelSagittaLogic_"+std::to_string(i),
                                                            0, 0, 0);
      m_DisplayAction->AddVolume(foamMidWheelSagittaLogic, "FoamEndWheel");

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

      new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, -layer_length_outer[i]/2 + foam_length / 2),
                        foamMidWheelSagittaLogic,
                        "foamMidWheelSagittaPlacedFrontTop"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, layer_length_outer[i]/2 - foam_length / 2),
                        foamMidWheelSagittaLogic,
                        "foamMidWheelSagittaPlacedBackTop"+std::to_string(i),
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

      new G4PVPlacement(G4Transform3D(bstlayer_rotm180, G4ThreeVector( 0.0, 0.0, -layer_length_outer[i]/2 + foam_length / 2)),
                        foamMidWheelSagittaLogic,
                        "foamMidWheelSagittaPlacedFrontBottom"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());
      new G4PVPlacement(G4Transform3D(bstlayer_rotm180, G4ThreeVector( 0.0, 0.0, layer_length_outer[i]/2 - foam_length / 2)),
                        foamMidWheelSagittaLogic,
                        "foamMidWheelSagittaPlacedBackBottom"+std::to_string(i),
                        mother,
                        0, 0, OverlapCheck());
    }
  }
  if(do_external_supports){

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
  }
  return;
}

G4Material* PHG4BSTDetector::MakeCarbonFoamMaterial(){
  G4Material* carbon_foam = GetDetectorMaterial("C_FOAM_BST", false);  // false suppresses warning that material does not exist
  if(!carbon_foam){
    G4double density;
    G4int ncomponents;
    carbon_foam = new G4Material("C_FOAM_BST", density = 0.50 * g / cm3, ncomponents = 2); // VERY CONSERVATIVE DENSITY
    // carbon_foam = new G4Material("C_FOAM_BST", density = 0.26 * g / cm3, ncomponents = 2); // CONSERVATIVE DENSITY
    // carbon_foam = new G4Material("C_FOAM_BST", density = 0.06 * g / cm3, ncomponents = 2); // LIGHTEST DENSITY
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
