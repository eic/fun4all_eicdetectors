#include "PHG4CylinderStripDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Element.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TMath.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>  // for operator<<, endl, basic_ost...
#include <map>
#include <sstream>

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________
PHG4CylinderStripDetector::PHG4CylinderStripDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_CylinderPhysicalVolume(nullptr)
  , m_Layer(lyr)
  , thickness(0)
{
  barwidth = 3 * mm;
}

//_______________________________________________________________
//bool PHG4CylinderStripDetector::IsInDetector(G4VPhysicalVolume *volume) const
//{
//  set<G4VPhysicalVolume *>::const_iterator iter =
//      m_PhysicalVolumesSet.find(volume);
//  if (iter != m_PhysicalVolumesSet.end())
//  {
//    return true;
//  }
//
//  return false;
//}

//_______________________________________________________________

bool PHG4CylinderStripDetector::IsInTileC(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter =
      m_CylinderCPhysicalVolume.find(volume);
  if (iter != m_CylinderCPhysicalVolume.end())
  {
    return true;
  }

  return false;
}

//_______________________________________________________________

bool PHG4CylinderStripDetector::IsInTileZ(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter =
      m_CylinderZPhysicalVolume.find(volume);
  if (iter != m_CylinderZPhysicalVolume.end())
  {
    return true;
  }

  return false;
}

//_______________________________________________________________
void PHG4CylinderStripDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  BuildMaterials();
  G4Material *TrackerMaterial = G4Material::GetMaterial(m_Params->get_string_param("gas"));

  if (!TrackerMaterial)
  {
    std::cout << "Error: Can not set Micromegas Gas" << std::endl;
    gSystem->Exit(1);
  }

  // components
  enum components
  {
    coverlay,
    CuGround,
    PCB,
    CuStrips,
    KaptonStrips,
    ResistiveStrips,
    Gas1,
    Mesh,
    Gas2,
    DriftCuElectrode,
    DriftKapton,
    DriftCuGround,
    kNcomponents
  };
  // component names
  G4String names[] = {
      "coverlay",
      "CuGround",
      "PCB",
      "CuStrips",
      "KaptonStrips",
      "ResistiveStrips",
      "Gas1",
      "Mesh",
      "Gas2",
      "DriftCuElectrode",
      "DriftKapton",
      "DriftCuGround"};

  // thicknesses
  map<int, float> thick;
  thick[coverlay] = 0.0050000 * cm;
  thick[CuGround] = 0.0001580 * cm;
  thick[PCB] = 0.0100000 * cm;
  thick[CuStrips] = 0.0012000 * cm;
  thick[KaptonStrips] = 0.0075000 * cm;
  thick[ResistiveStrips] = 0.0020000 * cm;
  //thick[ Gas1            ] = 0.0020000 * cm;
  thick[Gas1] = m_Params->get_double_param("gas1thickness") * cm;
  thick[Mesh] = 0.0018000 * cm;
  //thick[ Gas2            ] = 0.3000000 * cm;
  thick[Gas2] = m_Params->get_double_param("gas2thickness") * cm;
  thick[DriftCuElectrode] = 0.0005000 * cm;
  thick[DriftKapton] = 0.0250000 * cm;
  thick[DriftCuGround] = 0.0000410 * cm;

  // media
  map<int, G4Material *> media;
  media[coverlay] = G4Material::GetMaterial("myKapton");
  media[CuGround] = G4Material::GetMaterial("myCopper");
  media[PCB] = G4Material::GetMaterial("myFR4");
  media[CuStrips] = G4Material::GetMaterial("myMMStrips");
  media[KaptonStrips] = G4Material::GetMaterial("myKapton");
  media[ResistiveStrips] = G4Material::GetMaterial("myMMResistivePaste");
  media[Gas1] = TrackerMaterial;
  media[Mesh] = G4Material::GetMaterial("myMMMesh");
  media[Gas2] = TrackerMaterial;
  media[DriftCuElectrode] = G4Material::GetMaterial("myCopper");
  media[DriftKapton] = G4Material::GetMaterial("myKapton");
  media[DriftCuGround] = G4Material::GetMaterial("myCopper");

  // color
  map<int, G4Colour> color;
  color[coverlay] = G4Colour(204 / 255., 153 / 255., 0);
  color[CuGround] = G4Colour::Brown();
  color[PCB] = G4Colour::Green();
  color[CuStrips] = G4Colour::Brown();
  color[KaptonStrips] = G4Colour::Brown();
  color[ResistiveStrips] = G4Colour::Black();
  color[Gas1] = G4Colour::Grey();
  color[Mesh] = G4Colour::White();
  color[Gas2] = G4Colour::Grey();
  color[DriftCuElectrode] = G4Colour::Brown();
  color[DriftKapton] = G4Colour::Brown();
  color[DriftCuGround] = G4Colour(51 / 255., 26 / 255., 0);

  // components meca PCB
  enum components_meca
  {
    Cu1_meca,
    FR4_1_meca,
    Cu2_meca,
    FR4_2_meca,
    Cu3_meca,
    kNcomponents_meca
  };
  // component names
  G4String names_meca[] = {
      "Cu1_meca",
      "FR4_1_meca",
      "Cu2_meca",
      "FR4_2_meca",
      "Cu3_meca"};

  // thicknesses
  map<int, float> thick_meca;
  thick_meca[Cu1_meca] = 25 * um;
  thick_meca[FR4_1_meca] = 100 * um;
  thick_meca[Cu2_meca] = 25 * um;
  thick_meca[FR4_2_meca] = 100 * um;
  thick_meca[Cu3_meca] = 9 * um;

  // media
  map<int, G4Material *> media_meca;
  media_meca[Cu1_meca] = G4Material::GetMaterial("myCopper");
  media_meca[FR4_1_meca] = G4Material::GetMaterial("myFR4");
  media_meca[Cu2_meca] = G4Material::GetMaterial("myCopper");
  media_meca[FR4_2_meca] = G4Material::GetMaterial("myFR4");
  media_meca[Cu3_meca] = G4Material::GetMaterial("myCopper");

  // color
  map<int, G4Colour> color_meca;
  color_meca[Cu1_meca] = G4Colour::Brown();
  color_meca[FR4_1_meca] = G4Colour::Green();
  color_meca[Cu2_meca] = G4Colour::Brown();
  color_meca[FR4_2_meca] = G4Colour::Green();
  color_meca[Cu3_meca] = G4Colour::Brown();
  // determine length of cylinder using PHENIX's rapidity coverage if flag is true
  double radius = m_Params->get_double_param("radius") * cm;
  double gap = m_Params->get_double_param("gap") * cm;
  double gas_deadzone = m_Params->get_double_param("deadzone") * cm;
  int nhit = 1;

  int nCZlayer = 2;
  if (m_Params->get_int_param("use_2Dreadout"))
  {
    nhit = m_Params->get_int_param("nhit");
    nCZlayer = 1;
    gap = 0;
  }

  if (nCZlayer == 2)
    thick[CuStrips] = 9 * um;
  else if (nCZlayer == 1)
    thick[CuStrips] = 25 * um;

  thickness = 0;
  for (map<int, float>::iterator iter = thick.begin(); iter != thick.end(); ++iter)
  {
    thickness += iter->second;
  }
  // cout << "2D readout? " << m_Params->get_int_param("use_2Dreadout") << "  Cu strips thickness : " << thick[CuStrips] <<  endl;
  cout << "The tile thickness is " << thickness / mm << " mm" << endl;
  if (barwidth < thickness) barwidth = thickness;

  // make a "sub-world": a tube that contains all the tiles of this layer
  // --------------------------------------------------------------------

  // max thickness (cm)

  // radius is defined from the middle of the tracker layer, and now move it to the lower boundary of the bottom MM layer (1D), or to the lower boundary of the MM layer (2D)
  radius = radius - gap / 2. - thickness / 2. * nCZlayer;

  double spaceforhollowbar = 6 * mm;

  G4VSolid *cylinder_solid = new G4Tubs(G4String(GetName()),
                                        radius - 0.001 * mm - spaceforhollowbar / 2,
                                        radius + thickness * nCZlayer + gap + 0.001 * mm + spaceforhollowbar / 2,
                                        m_Params->get_double_param("length") * cm / 2. + barwidth + 2 * mm, 0, twopi);
  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid,
                                                        G4Material::GetMaterial("myAir"),
                                                        G4String(GetName()));
  G4VisAttributes *vis = new G4VisAttributes(G4Color(G4Colour::Grey()));  // grey is good to see the tracks in the display
  vis->SetForceSolid(false);
  vis->SetVisibility(false);
  vis->SetForceLineSegmentsPerCircle(30);
  cylinder_logic->SetVisAttributes(vis);

  // add mother logical volume to subsystem
  PHG4Subsystem *mysys = GetMySubsystem();
  mysys->SetLogicalVolume(cylinder_logic);

  // compute how many tiles
  // ----------------------

  float maxTileWidth = 50. * cm;  // cm
  float spacer = 2. * cm;         // cm, mechanical space between tiles, to be filled with carbon fiber

  float circumference = twopi * radius;
  int Ntiles = ceil(circumference / (maxTileWidth + spacer));

  float tileW = (circumference - Ntiles * spacer) / Ntiles;

  cout << setw(10) << radius << setw(10) << "Ntiles: " << Ntiles << setw(20) << "tileW: " << tileW << endl;

  // make 1 tile
  // -----------

  float deltaPhi = tileW / radius * TMath::RadToDeg();
  spacer = radius * 360. / Ntiles * deg - tileW;
  double steplimits = m_Params->get_double_param("steplimits") * cm;
  G4UserLimits *g4userlimits = nullptr;
  if (isfinite(steplimits))
  {
    g4userlimits = new G4UserLimits(steplimits);
  }

  G4VSolid *tile_o = new G4Tubs(G4String(GetName()) + "_tile",
                                radius - 0.001 * mm - spaceforhollowbar / 2,
                                radius + thickness * nCZlayer + gap + 0.001 * mm + spaceforhollowbar / 2,
                                m_Params->get_double_param("length") * cm / 2. + barwidth + 2 * mm, 0, 360. / Ntiles * deg - 0.01 * deg);
  G4LogicalVolume *tile_o_logic = new G4LogicalVolume(tile_o,
                                                      G4Material::GetMaterial("myAir"),
                                                      G4String(GetName()) + "_tile_logic");
  tile_o_logic->SetVisAttributes(vis);

  double thickness_mecaPCB = 25 * um + 100 * um + 25 * um + 100 * um + 9 * um;
  double radius_mecaPCB = radius + thickness / 2. - thickness_mecaPCB / 2.;
  // build a mother volume filled by carbon fiber using union volume
  G4VSolid *mecaPCB_solid = new G4Tubs("mecaPCB_solid",
                                       radius_mecaPCB - 0.1 * um,
                                       radius_mecaPCB + thickness_mecaPCB + 0.1 * um,
                                       m_Params->get_double_param("length") * cm / 2.,
                                       360. / Ntiles * deg - (spacer - 2 * barwidth + 1 * mm) / radius_mecaPCB * radian,
                                       (spacer - 2 * barwidth) / (radius_mecaPCB) *radian);
  G4VSolid *MM_solid = new G4Tubs("MM_solid",
                                  radius - 0.1 * um,
                                  radius + thickness + 0.1 * um,
                                  m_Params->get_double_param("length") * cm / 2. + barwidth + .1 * mm,
                                  barwidth / radius * radian - 0.5 * mm / radius * radian,
                                  deltaPhi * deg + 0.5 * mm / radius * radian);
  G4VSolid *bar_solid = GetHollowBar();

  G4RotationMatrix *zrot = new G4RotationMatrix;
  double rotang_bar2 = barwidth / 2. / radius * radian;
  zrot->rotateZ(rotang_bar2);
  G4VSolid *u1 = new G4UnionSolid("MM+bar1", MM_solid, bar_solid, G4Transform3D(*zrot, G4ThreeVector((radius + thickness / 2.) * TMath::Cos(rotang_bar2), (radius + thickness / 2.) * TMath::Sin(rotang_bar2), 0)));

  rotang_bar2 = 360. / Ntiles * deg + rotang_bar2 - ((spacer - barwidth) / radius * radian) - 0.5 * mm / radius * radian;
  zrot->setDelta(twopi - (360. / Ntiles * deg - ((spacer - barwidth) / radius * radian) - 90 * deg));
  u1 = new G4UnionSolid("MM+bar1+bar2", u1, bar_solid, zrot, G4ThreeVector((radius + thickness / 2.) * TMath::Cos(rotang_bar2), (radius + thickness / 2.) * TMath::Sin(rotang_bar2), 0));
  //G4RotationMatrix* zrotx = new G4RotationMatrix();
  //zrotx->rotateZ(deltaPhi*deg + 2*barwidth/radius *radian + 2*mm/radius);
  //u1 = new G4UnionSolid("MM+bar1+bar2+arch1+arch2+meca", u1, mecaPCB_solid, G4Transform3D(*zrotx, G4ThreeVector(0,0,0)));
  u1 = new G4UnionSolid("MM+bar1+bar2+arch1+arch2+meca", u1, mecaPCB_solid);
  // logic volume filled with carbon fiber
  G4LogicalVolume *u1_C_logic = new G4LogicalVolume(u1,
                                                    G4Material::GetMaterial("myCfiber"),
                                                    "u1_C_logic");

  vis = new G4VisAttributes(G4Color(G4Colour::Grey()));  // grey is good to see the tracks in the display
  vis->SetForceSolid(false);
  vis->SetVisibility(true);
  vis->SetForceLineSegmentsPerCircle(30);
  u1_C_logic->SetVisAttributes(vis);
  zrot->setDelta(0);
  new G4PVPlacement(G4Transform3D(*zrot, G4ThreeVector(0, 0, 0)),
                    u1_C_logic,
                    "u1_C_phys",
                    tile_o_logic, false, 0, OverlapCheck());

  float Rm = radius;
  float RM = Rm;
  float Rm_meca = radius_mecaPCB;
  float RM_meca = Rm_meca;

  G4VSolid *tile_o_comp = nullptr;
  G4LogicalVolume *tile_o_comp_logic = nullptr;

  for (int ic = 0; ic < kNcomponents; ic++)
  {
    G4UserLimits *g4userlimits_gas = nullptr;
    if (ic == Gas2)
      g4userlimits_gas = g4userlimits;

    G4String cname = G4String(GetName()) + "_tileC" + "_" + names[ic];

    zrot->setDelta(barwidth / radius * radian);
    G4VPhysicalVolume *phys = nullptr;

    vis = new G4VisAttributes(G4Color(color[ic]));  // grey is good to see the tracks in the display
    vis->SetForceSolid(true);
    vis->SetVisibility(true);
    vis->SetForceLineSegmentsPerCircle(30);

    if (ic == Gas2)
    {
      thick[ic] /= nhit;
      for (int ii = 0; ii < nhit; ii++)
      {
        RM = Rm + thick[ic];
        tile_o_comp = new G4Tubs(cname + "_solid",
                                 Rm,
                                 RM,
                                 m_Params->get_double_param("length") * cm / 2. - gas_deadzone, gas_deadzone / Rm * radian, deltaPhi * deg - 2 * gas_deadzone / Rm * radian);
        tile_o_comp_logic = new G4LogicalVolume(tile_o_comp,
                                                media[ic],
                                                cname + "_logic",
                                                nullptr,
                                                nullptr,
                                                g4userlimits_gas);
        tile_o_comp_logic->SetVisAttributes(vis);
        Rm = RM;
        phys = new G4PVPlacement(G4Transform3D(*zrot, G4ThreeVector(0, 0, 0)),
                                 tile_o_comp_logic,
                                 cname + "_phys",
                                 u1_C_logic, false, 0, OverlapCheck());
      }
    }
    else
    {
      RM = Rm + thick[ic];
      tile_o_comp = new G4Tubs(cname + "_solid",
                               Rm,
                               RM,
                               m_Params->get_double_param("length") * cm / 2., 0, deltaPhi * deg);
      tile_o_comp_logic = new G4LogicalVolume(tile_o_comp,
                                              media[ic],
                                              cname + "_logic",
                                              nullptr,
                                              nullptr,
                                              g4userlimits_gas);
      tile_o_comp_logic->SetVisAttributes(vis);
      Rm = RM;
      phys = new G4PVPlacement(G4Transform3D(*zrot, G4ThreeVector(0, 0, 0)),
                               tile_o_comp_logic,
                               cname + "_phys",
                               u1_C_logic, false, 0, OverlapCheck());
    }
    if (ic == Gas2)
      m_CylinderCPhysicalVolume.insert(phys);
  }

  for (int ic = 0; ic < kNcomponents_meca; ic++)
  {
    G4String cname = G4String(GetName()) + "_tileC" + "_" + names_meca[ic];

    RM_meca = Rm_meca + thick_meca[ic];

    tile_o_comp = new G4Tubs(cname + "_solid",
                             Rm_meca,
                             RM_meca,
                             m_Params->get_double_param("length") * cm / 2.,
                             0,
                             (spacer - 2 * barwidth - 0.001 * mm) / (radius_mecaPCB) *radian);
    tile_o_comp_logic = new G4LogicalVolume(tile_o_comp,
                                            media_meca[ic],
                                            cname + "_logic");
    vis = new G4VisAttributes(G4Color(color_meca[ic]));  // grey is good to see the tracks in the display
    vis->SetForceSolid(true);
    vis->SetVisibility(true);
    vis->SetForceLineSegmentsPerCircle(30);
    tile_o_comp_logic->SetVisAttributes(vis);
    G4RotationMatrix *zrot_tmp = new G4RotationMatrix();
    zrot_tmp->rotateZ(360. / Ntiles * deg - (spacer - 2 * barwidth + 1 * mm) / radius_mecaPCB * radian);
    new G4PVPlacement(G4Transform3D(*zrot_tmp,
                                    G4ThreeVector(0, 0, 0)),
                      tile_o_comp_logic,
                      cname + "_phys",
                      u1_C_logic,
                      false,
                      0,
                      OverlapCheck());
    Rm_meca = RM_meca;
  }

  if (nCZlayer == 2)
  {
    //Rm += gap;
    //RM += gap;
    //Rm_meca = Rm + thickness/2. - thickness_mecaPCB/2.;
    //RM_meca = Rm_meca;
    Rm = radius;
    RM = Rm;
    Rm_meca = radius_mecaPCB;
    RM_meca = Rm_meca;
    // logic volume filled with carbon fiber
    G4LogicalVolume *u1_Z_logic = new G4LogicalVolume(u1,
                                                      G4Material::GetMaterial("myCfiber"),
                                                      "u1_Z_logic");
    vis = new G4VisAttributes(G4Color(G4Colour::Grey()));  // grey is good to see the tracks in the display
    vis->SetForceSolid(false);
    vis->SetVisibility(true);
    vis->SetForceLineSegmentsPerCircle(30);
    u1_Z_logic->SetVisAttributes(vis);
    new G4PVPlacement(0, G4ThreeVector((gap + thickness) * TMath::Cos(360. / Ntiles * deg / 2.), (gap + thickness) * TMath::Sin(360. / Ntiles * deg / 2.), 0),
                      u1_Z_logic,
                      "u1_Z_phys",
                      tile_o_logic, false, 0, OverlapCheck());
    for (int ic = 0; ic < kNcomponents; ic++)
    {
      G4UserLimits *g4userlimits_gas = nullptr;
      if (ic == Gas2)
        g4userlimits_gas = g4userlimits;

      G4String cname = G4String(GetName()) + "_tileZ" + "_" + names[ic];

      RM = Rm + thick[ic];

      if (ic == Gas2)
      {
        thick[ic] /= nhit;
        for (int ii = 0; ii < nhit; ii++)
        {
          RM = Rm + thick[ic];
          tile_o_comp = new G4Tubs(cname + "_solid",
                                   Rm,
                                   RM,
                                   m_Params->get_double_param("length") * cm / 2. - gas_deadzone, gas_deadzone / Rm * radian, deltaPhi * deg - 2 * gas_deadzone / Rm * radian);
          tile_o_comp_logic = new G4LogicalVolume(tile_o_comp,
                                                  media[ic],
                                                  cname + "_logic",
                                                  nullptr,
                                                  nullptr,
                                                  g4userlimits_gas);
          Rm = RM;
        }
      }
      else
      {
        RM = Rm + thick[ic];
        tile_o_comp = new G4Tubs(cname + "_solid",
                                 Rm,
                                 RM,
                                 m_Params->get_double_param("length") * cm / 2., 0, deltaPhi * deg);
        tile_o_comp_logic = new G4LogicalVolume(tile_o_comp,
                                                media[ic],
                                                cname + "_logic",
                                                nullptr,
                                                nullptr,
                                                g4userlimits_gas);
        Rm = RM;
      }
      vis = new G4VisAttributes(G4Color(color[ic]));  // grey is good to see the tracks in the display
      vis->SetForceSolid(true);
      vis->SetVisibility(true);
      vis->SetForceLineSegmentsPerCircle(30);
      tile_o_comp_logic->SetVisAttributes(vis);
      zrot->setDelta(barwidth / radius * radian);
      G4VPhysicalVolume *phys = new G4PVPlacement(G4Transform3D(*zrot, G4ThreeVector(0, 0, 0)),
                                                  tile_o_comp_logic,
                                                  cname + "_phys",
                                                  u1_Z_logic, false, 0, OverlapCheck());
      if (ic == Gas2)
        m_CylinderZPhysicalVolume.insert(phys);
    }
    for (int ic = 0; ic < kNcomponents_meca; ic++)
    {
      G4String cname = G4String(GetName()) + "_tileZ" + "_" + names_meca[ic];

      RM_meca = Rm_meca + thick_meca[ic];

      tile_o_comp = new G4Tubs(cname + "_solid",
                               Rm_meca,
                               RM_meca,
                               m_Params->get_double_param("length") * cm / 2.,
                               0,
                               (spacer - 2 * barwidth - 0.001 * mm) / (radius_mecaPCB) *radian);
      tile_o_comp_logic = new G4LogicalVolume(tile_o_comp,
                                              media_meca[ic],
                                              cname + "_logic");
      vis = new G4VisAttributes(G4Color(color_meca[ic]));  // grey is good to see the tracks in the display
      vis->SetForceSolid(true);
      vis->SetVisibility(true);
      vis->SetForceLineSegmentsPerCircle(30);
      tile_o_comp_logic->SetVisAttributes(vis);
      G4RotationMatrix *zrot_tmp2 = new G4RotationMatrix();
      zrot_tmp2->rotateZ(360. / Ntiles * deg - (spacer - 2 * barwidth + 1 * mm) / radius_mecaPCB * radian);
      new G4PVPlacement(G4Transform3D(*zrot_tmp2,
                                      G4ThreeVector(0, 0, 0)),
                        tile_o_comp_logic,
                        cname + "_phys",
                        u1_Z_logic,
                        false,
                        0,
                        OverlapCheck());
      Rm_meca = RM_meca;
    }
  }

  // repeate N tiles
  for (int i = 0; i < Ntiles; i++)
  {
    G4RotationMatrix *yRot = new G4RotationMatrix;
    yRot->rotateZ(i * 360. / Ntiles * deg);
    new G4PVPlacement(G4Transform3D(*yRot, G4ThreeVector(0, 0, 0)),
                      tile_o_logic,
                      G4String(GetName()) + "_tile_phys",
                      cylinder_logic,
                      true,
                      i,
                      OverlapCheck());
  }

  double phi0 = m_Params->get_double_param("phi0") * deg;
  G4RotationMatrix *yRot0 = new G4RotationMatrix;
  yRot0->rotateZ(phi0);
  m_CylinderPhysicalVolume = new G4PVPlacement(G4Transform3D(*yRot0, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm)),
                                               cylinder_logic,
                                               G4String(GetName()),
                                               logicWorld, 0, false, OverlapCheck());
}

void PHG4CylinderStripDetector::BuildMaterials()
{
  // get the list of chemical elements
  // ---------------------------------

  //G4Element *N  = G4Element::GetElement( "N" );
  //G4Element *O  = G4Element::GetElement( "O" );
  //G4Element *H  = G4Element::GetElement( "H" );
  //G4Element *C  = G4Element::GetElement( "C" );

  // get the list of NIST materials
  // ---------------------------------
  G4Material *G4_N = G4Material::GetMaterial("G4_N");
  G4Material *G4_O = G4Material::GetMaterial("G4_O");
  G4Material *G4_C = G4Material::GetMaterial("G4_C");
  G4Material *G4_H = G4Material::GetMaterial("G4_H");
  G4Material *G4_Si = G4Material::GetMaterial("G4_Si");
  G4Material *G4_Ar = G4Material::GetMaterial("G4_Ar");
  G4Material *G4_Cr = G4Material::GetMaterial("G4_Cr");
  G4Material *G4_Fe = G4Material::GetMaterial("G4_Fe");
  G4Material *G4_Mn = G4Material::GetMaterial("G4_Mn");
  G4Material *G4_Ni = G4Material::GetMaterial("G4_Ni");
  G4Material *G4_Cu = G4Material::GetMaterial("G4_Cu");

  // combine elements
  // ----------------
  G4int ncomponents;
  G4double fraction;
  G4double temperature = 298.15 * kelvin;
  G4double pressure = 1. * atmosphere;

  // air
  if (!G4Material::GetMaterial("myAir", false))
  {
    G4Material *myAir = new G4Material("myAir", 0.001205 * g / cm3, ncomponents = 2, kStateGas, temperature, pressure);
    myAir->AddMaterial(G4_N, fraction = 0.77);
    myAir->AddMaterial(G4_O, fraction = 0.23);
  }

  // FR4
  if (!G4Material::GetMaterial("myFR4", false))
  {
    G4Material *myFR4 = new G4Material("myFR4", 1.860 * g / cm3, ncomponents = 4, kStateSolid);
    myFR4->AddMaterial(G4_C, fraction = 0.43550);
    myFR4->AddMaterial(G4_H, fraction = 0.03650);
    myFR4->AddMaterial(G4_O, fraction = 0.28120);
    myFR4->AddMaterial(G4_Si, fraction = 0.24680);
  }

  // Kapton
  if (!G4Material::GetMaterial("myKapton", false))
  {
    G4Material *myKapton = new G4Material("myKapton", 1.420 * g / cm3, ncomponents = 4, kStateSolid);
    myKapton->AddMaterial(G4_C, 0.6911330);
    myKapton->AddMaterial(G4_H, 0.0263620);
    myKapton->AddMaterial(G4_N, 0.0732700);
    myKapton->AddMaterial(G4_O, 0.2092350);
  }

  // MMgas
  if (!G4Material::GetMaterial("myMMGas", false))
  {
    G4Material *myMMGas = new G4Material("myMMGas", 0.00170335 * g / cm3, ncomponents = 3, kStateGas, temperature, pressure);
    myMMGas->AddMaterial(G4_Ar, 0.900);
    myMMGas->AddMaterial(G4_C, 0.0826586);
    myMMGas->AddMaterial(G4_H, 0.0173414);
  }

  // MMMesh
  if (!G4Material::GetMaterial("myMMMesh", false))
  {
    G4Material *myMMMesh = new G4Material("myMMMesh", 2.8548 * g / cm3, ncomponents = 5, kStateSolid);
    myMMMesh->AddMaterial(G4_Cr, 0.1900);
    myMMMesh->AddMaterial(G4_Fe, 0.6800);
    myMMMesh->AddMaterial(G4_Mn, 0.0200);
    myMMMesh->AddMaterial(G4_Ni, 0.1000);
    myMMMesh->AddMaterial(G4_Si, 0.0100);
  }

  // MMStrips
  if (!G4Material::GetMaterial("myMMStrips", false))
  {
    G4Material *myMMStrips = new G4Material("myMMStrips", 5.248414 * g / cm3, G4_Cu, kStateSolid);
    cout << myMMStrips->GetName() << endl;
  }

  // MMResistivePaste
  if (!G4Material::GetMaterial("myMMResistivePaste", false))
  {
    G4Material *myMMResistivePaste = new G4Material("myMMResistivePaste", 0.77906 * g / cm3, G4_C, kStateSolid);
    cout << myMMResistivePaste->GetName() << endl;
  }

  // Copper
  if (!G4Material::GetMaterial("myCopper", false))
  {
    G4Material *myCopper = new G4Material("myCopper", 8.9600 * g / cm3, G4_Cu, kStateSolid);
    cout << myCopper->GetName() << endl;
  }

  // Carbon fiber
  if (!G4Material::GetMaterial("myCfiber", false))
  {
    G4Material *myCfiber = new G4Material("myCfiber", 1.80 * g / cm3, G4_C, kStateSolid);
    cout << myCfiber->GetName() << endl;
  }
}

G4VSolid *PHG4CylinderStripDetector::GetHollowBar()
{
  double length = m_Params->get_double_param("length") * cm;
  G4VSolid *box1 = new G4Box("bar0", barwidth / 2, barwidth / 2, (length + barwidth * 2) / 2);
  G4VSolid *box2 = new G4Box("bar1", 1 * mm / 2, 1 * mm / 2, (length + barwidth * 2) / 2 + 1 * mm);
  G4VSolid *bar = new G4SubtractionSolid("HollowBar", box1, box2);
  return bar;
}
