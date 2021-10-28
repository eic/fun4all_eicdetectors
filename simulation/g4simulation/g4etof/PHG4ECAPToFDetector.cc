#include "PHG4ECAPToFDetector.h"
#include <g4main/PHG4Detector.h>  // for PHG4Detector
#include <phparameter/PHParameters.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4Element.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVDivision.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>  // for operator<<, std::endl, bas...
#include <sstream>
#include <string>

class PHCompositeNode;
class PHG4Subsystem;

PHG4ECAPToFDetector::PHG4ECAPToFDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , Phys(nullptr)
  , fhc_phys(nullptr)
  , fpcb_phys(nullptr)
  , fpcbcu_phys(nullptr)
  , fmylar_phys(nullptr)
  , fcarbon_phys(nullptr)
  , fglass_phys()
  , fgas_phys()
  , mcarbon_phys(nullptr)
  , mmylar_phys(nullptr)
  , mpcbcu_phys(nullptr)
  , mpcb_phys(nullptr)
  , mpcbcu2_phys(nullptr)
  , mmylar2_phys(nullptr)
  , mcarbon2_phys(nullptr)
  , bglass_phys()
  , bgas_phys()
  , bhc_phys(nullptr)
  , bpcb_phys(nullptr)
  , bpcbcu_phys(nullptr)
  , bmylar_phys(nullptr)
  , bcarbon_phys(nullptr)
  , m_Active(m_Params->get_int_param("active"))
  , m_Layer(lyr)
{
}

int PHG4ECAPToFDetector::IsInToF(const G4VPhysicalVolume *volume) const
{
  if (m_Active)
  {
    for (int i = 0; i < 6; i++)  // neeed to define # of gas layers
    {
      if (volume == fgas_phys[i] || volume == bgas_phys[i])
      {
        return 1;
      }
    }
  }
  return 0;
}

void PHG4ECAPToFDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  //gaps and thickness
  int n_fgas_layer = m_Params->get_int_param("n_fgas_layer");
  int n_bgas_layer = m_Params->get_int_param("n_bgas_layer");
  double gas_gap = m_Params->get_double_param("gas_gap") * cm;
  double glass_thick = m_Params->get_double_param("glass_thick") * cm;
  double carbon_thick = m_Params->get_double_param("Carbon_thick") * cm;
  double pcb_thick = m_Params->get_double_param("pcb_thick") * cm;
  double cu_thick = m_Params->get_double_param("cu_thick") * cm;
  double honeycomb_thick = m_Params->get_double_param("honeycomb_thick") * cm;
  double mylar_thick = m_Params->get_double_param("mylar_thick") * cm;
  double Rin = m_Params->get_double_param("Rin") * cm;
  double Rout = m_Params->get_double_param("Rout") * cm;
  double z_begin = m_Params->get_double_param("z_begin") * cm;

  if (Verbosity() > 1)
  {
    std::cout << " passed on parameters from macros :: " << std::endl;
    std::cout << "n_fgas_layer :" << n_fgas_layer << std::endl;
    std::cout << "n_bgas_layer : " << n_bgas_layer << std::endl;
    std::cout << "gas_gap : " << gas_gap << std::endl;
    std::cout << "glass_thick : " << glass_thick << std::endl;
    std::cout << "carbon_thick : " << carbon_thick << std::endl;
    std::cout << "pcb_thick : " << pcb_thick << std::endl;
    std::cout << "cu_thick : " << cu_thick << std::endl;
    std::cout << "honeycomb_thick : " << honeycomb_thick << std::endl;
    std::cout << "mylar_thick : " << mylar_thick << std::endl;
    std::cout << "Rin : " << Rin << std::endl;
    std::cout << "Rout : " << Rout << std::endl;
    std::cout << "z_begin : " << z_begin << std::endl;
  }

  double tot_thick = (n_fgas_layer + n_bgas_layer) * gas_gap + (n_fgas_layer + n_bgas_layer + 2) * glass_thick + 4. * carbon_thick + 3. * pcb_thick + 4. * cu_thick + 4 * mylar_thick + 2. * honeycomb_thick;
  double posz = z_begin + 0.5 * tot_thick;

  if (Verbosity() > 1)
  {
    std::cout << "Mother vol dimensions : " << std::endl;
    std::cout << " z begin :" << z_begin << std::endl;
    std::cout << " tot_thick :" << tot_thick << std::endl;
    std::cout << " mid pos :" << posz << std::endl;
  }
  G4Material *tof_mother_mat = G4Material::GetMaterial(m_Params->get_string_param("material"));

  G4Tubs *Solid = new G4Tubs("ToF_GVol_Solid", Rin, Rout, tot_thick / 2., 0., 360. * deg);
  G4LogicalVolume *Logic = new G4LogicalVolume(Solid, tof_mother_mat, "ToF_GVol_Logic", 0, 0, 0);
  G4VisAttributes *attr_ToF_GVol = new G4VisAttributes(G4Color(0.3, 0.5, 0.9, 0.9));
  //attr_ToF_GVol->SetColor(G4Color::Green());
  attr_ToF_GVol->SetForceSolid(true);
  Logic->SetVisAttributes(attr_ToF_GVol);
  Phys = new G4PVPlacement(0, G4ThreeVector(0, 0, posz), Logic, "E_CAP_ToF_Physics", logicWorld, 0, false, OverlapCheck());

  G4Material *G4_O = G4Material::GetMaterial("G4_O");
  G4Material *G4_Si = G4Material::GetMaterial("G4_Si");
  G4Material *G4_Na = G4Material::GetMaterial("G4_Na");
  G4Material *G4_Ca = G4Material::GetMaterial("G4_Ca");
  G4Material *G4_Cu = G4Material::GetMaterial("G4_Cu");
  G4Material *G4_FR4 = G4Material::GetMaterial("FR4");
  G4Material *G4_Nomex = G4Material::GetMaterial("NOMEX");
  G4Material *G4_Mylar = G4Material::GetMaterial("G4_MYLAR");
  G4Material *G4_gas = G4Material::GetMaterial("G4_Ar");

  //Plate Glass material
  G4Material *plateglass = new G4Material("plateglass", 2.4 * g / cm3, 4, kStateSolid);
  plateglass->AddMaterial(G4_O, 0.4598);
  plateglass->AddMaterial(G4_Na, 0.096441);
  plateglass->AddMaterial(G4_Si, 0.336553);
  plateglass->AddMaterial(G4_Ca, 0.107205);

  //Now define various elements from front of ToF

  //Front Honeycomb
  G4Tubs *fhc_solid = new G4Tubs("front_honeycomb_solid", Rin, Rout, 0.5 * honeycomb_thick, 0., 360 * deg);
  G4LogicalVolume *fhc_logic = new G4LogicalVolume(fhc_solid, G4_Nomex, "front_honeycomb_Logic", 0, 0, 0);

  G4VisAttributes *tof_hc_att = new G4VisAttributes();
  tof_hc_att->SetColour(G4Colour::Yellow());
  tof_hc_att->SetVisibility(true);
  tof_hc_att->SetForceSolid(true);
  fhc_logic->SetVisAttributes(tof_hc_att);
  double fhc_pos = -0.5 * tot_thick + 0.5 * honeycomb_thick;

  fhc_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, fhc_pos), fhc_logic, "front_hc_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Front HC mid pos :" << fhc_pos << std::endl;
  //Front PCB

  G4Tubs *fpcb_solid = new G4Tubs("front_pcb_solid", Rin, Rout, 0.5 * pcb_thick, 0., 360 * deg);
  G4LogicalVolume *fpcb_logic = new G4LogicalVolume(fpcb_solid, G4_FR4, "front_PCB_Logic", 0, 0, 0);

  G4VisAttributes *tof_pcb_att = new G4VisAttributes();
  tof_pcb_att->SetColour(G4Colour::Green());
  tof_pcb_att->SetVisibility(true);
  tof_pcb_att->SetForceSolid(true);
  fpcb_logic->SetVisAttributes(tof_pcb_att);
  double fpcb_pos = fhc_pos + 0.5 * honeycomb_thick + 0.5 * pcb_thick;

  fpcb_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, fpcb_pos), fpcb_logic, "front_pcb_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Front PCB mid pos : " << fpcb_pos << std::endl;

  //Front PCB Cu
  G4Tubs *fpcbcu_solid = new G4Tubs("front_pcbcu_solid", Rin, Rout, 0.5 * cu_thick, 0., 360 * deg);
  G4LogicalVolume *fpcbcu_logic = new G4LogicalVolume(fpcbcu_solid, G4_Cu, "front_PCBCu_Logic", 0, 0, 0);

  G4VisAttributes *cu_att = new G4VisAttributes();
  cu_att->SetColour(G4Colour::Yellow());
  cu_att->SetVisibility(true);
  cu_att->SetForceSolid(true);
  fpcbcu_logic->SetVisAttributes(cu_att);
  double fpcbcu_pos = fpcb_pos + 0.5 * pcb_thick + 0.5 * cu_thick;

  fpcbcu_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, fpcbcu_pos), fpcbcu_logic, "front_pcbcu_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Front PCB Cu pos : " << fpcbcu_pos << std::endl;

  //Front Mylar
  G4Tubs *fmylar_solid = new G4Tubs("front_mylar_solid", Rin, Rout, 0.5 * mylar_thick, 0., 360 * deg);
  G4LogicalVolume *fmylar_logic = new G4LogicalVolume(fmylar_solid, G4_Mylar, "front_Mylar_Logic", 0, 0, 0);

  G4VisAttributes *mylar_att = new G4VisAttributes();
  mylar_att->SetColour(G4Colour::White());
  mylar_att->SetVisibility(true);
  mylar_att->SetForceSolid(true);
  fmylar_logic->SetVisAttributes(mylar_att);
  double fmylar_pos = fpcbcu_pos + 0.5 * cu_thick + 0.5 * mylar_thick;

  fmylar_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, fmylar_pos), fmylar_logic, "front_mylar_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Front mylar pos :" << fmylar_pos << std::endl;
  //Carbon layer
  G4Tubs *fcarbon_solid = new G4Tubs("front_carbon_solid", Rin, Rout, 0.5 * carbon_thick, 0., 360 * deg);
  G4LogicalVolume *fcarbon_logic = new G4LogicalVolume(fcarbon_solid, G4_Ca, "front_Carbon_Logic", 0, 0, 0);

  G4VisAttributes *carbon_att = new G4VisAttributes();
  carbon_att->SetColour(G4Colour::Black());
  carbon_att->SetVisibility(true);
  carbon_att->SetForceSolid(true);
  fcarbon_logic->SetVisAttributes(carbon_att);
  double fcarbon_pos = fmylar_pos + 0.5 * mylar_thick + 0.5 * carbon_thick;

  fcarbon_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, fcarbon_pos), fcarbon_logic, "front_carbon_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " front Carbon pos :" << fcarbon_pos << std::endl;
  //mrpc and gas layers
  double glass_begin_pos = fcarbon_pos + 0.5 * carbon_thick;
  double glass_mid_pos = 0.;
  double gas_begin_pos = glass_begin_pos + glass_thick;
  double gas_mid_pos = 0.;

  G4VisAttributes *glass_att = new G4VisAttributes();
  glass_att->SetColor(G4Colour::Magenta());
  glass_att->SetVisibility(true);
  glass_att->SetForceSolid(true);

  G4VisAttributes *gas_att = new G4VisAttributes();
  gas_att->SetColor(G4Colour::Blue());
  gas_att->SetVisibility(true);
  gas_att->SetForceSolid(true);

  G4Tubs *fglass_solid[n_fgas_layer + 1];
  G4LogicalVolume *fglass_logic[n_fgas_layer + 1];
  G4Tubs *fgas_solid[n_fgas_layer];
  G4LogicalVolume *fgas_logic[n_fgas_layer];

  for (int l = 0; l < n_fgas_layer + 1; l++)
  {
    /*
      fglass_solid[l] = new G4Tubs(Form("front_glass_solid_%d", l), Rin, Rout, 0.5*glass_thick, 0, 360.*deg);
      fglass_logic[l] = new G4LogicalVolume(fglass_solid[l], plateglass, Form("front_glass_Logic_%d", l),0 , 0,0);
    
      fglass_solid[l] = new G4Tubs("front_glass_solid_%d", Rin, Rout, 0.5*glass_thick, 0, 360.*deg);
      fglass_logic[l] = new G4LogicalVolume(fglass_solid[l], plateglass,"front_glass_Logic_%d",0 , 0,0);
    */

    std::string gl_solid_name = "fglass_solid_" + std::to_string(l);
    std::string gl_logic_name = "fglass_logic_" + std::to_string(l);
    std::string gl_phys_name = "fglass_phys_" + std::to_string(l);

    fglass_solid[l] = new G4Tubs(gl_solid_name, Rin, Rout, 0.5 * glass_thick, 0, 360. * deg);
    fglass_logic[l] = new G4LogicalVolume(fglass_solid[l], plateglass, gl_logic_name, 0, 0, 0);
    fglass_logic[l]->SetVisAttributes(glass_att);
    glass_mid_pos = glass_begin_pos + (l + 0.5) * glass_thick + l * gas_gap;

    //fglass_phys[l] =  new G4PVPlacement(0, G4ThreeVector(0, 0, glass_mid_pos), fglass_logic[l], sprintf("MRPC_layer_phys_%d",l), Logic, 0, false, OverlapCheck());
    fglass_phys[l] = new G4PVPlacement(0, G4ThreeVector(0, 0, glass_mid_pos), fglass_logic[l], gl_phys_name, Logic, 0, false, OverlapCheck());

    if (l < n_fgas_layer)
    {
      /*
	fgas_solid[l] = new G4Tubs(sprintf("front_gas_solid_%d", l), Rin, Rout, 0.5*gas_gap, 0, 360.*deg);
	fgas_logic[l] = new G4LogicalVolume(fgas_solid[l], G4_gas, sprintf("front_gas_Logic_%d", l),0 , 0,0);
      
	fgas_solid[l] = new G4Tubs("front_gas_solid_%d", Rin, Rout, 0.5*gas_gap, 0, 360.*deg);
	fgas_logic[l] = new G4LogicalVolume(fgas_solid[l], G4_gas,"front_gas_Logic_%d",0 , 0,0);
      */

      std::string gas_solid_name = "fgas_solid_" + std::to_string(l);
      std::string gas_logic_name = "fgas_logic_" + std::to_string(l);
      std::string gas_phys_name = "fgas_phys_" + std::to_string(l);

      fgas_solid[l] = new G4Tubs(gas_solid_name, Rin, Rout, 0.5 * gas_gap, 0, 360. * deg);
      fgas_logic[l] = new G4LogicalVolume(fgas_solid[l], G4_gas, gas_logic_name, 0, 0, 0);
      fgas_logic[l]->SetVisAttributes(gas_att);
      gas_mid_pos = gas_begin_pos + (l + 0.5) * gas_gap + l * glass_thick;

      //fgas_phys[l] =  new G4PVPlacement(0, G4ThreeVector(0, 0, gas_mid_pos), fgas_logic[l], sprintf("Gas_layer_phys_%d",l), Logic, 0, false, OverlapCheck());
      fgas_phys[l] = new G4PVPlacement(0, G4ThreeVector(0, 0, gas_mid_pos), fgas_logic[l], gas_phys_name, Logic, 0, false, OverlapCheck());
    }
    if (Verbosity() > 1)
    {
      std::cout << " Front glass layer : " << l << " glass mid : " << glass_mid_pos << std::endl;
      std::cout << "Front gas layer : " << l << " gas mid : " << gas_mid_pos << std::endl;
    }
  }

  //Middle Carbon layer

  G4Tubs *mcarbon_solid = new G4Tubs("mid_carbon_solid", Rin, Rout, 0.5 * carbon_thick, 0., 360 * deg);
  G4LogicalVolume *mcarbon_logic = new G4LogicalVolume(mcarbon_solid, G4_Ca, "mid_Carbon_Logic", 0, 0, 0);
  mcarbon_logic->SetVisAttributes(carbon_att);
  double mcarbon_pos = glass_mid_pos + 0.5 * glass_thick + 0.5 * carbon_thick;
  mcarbon_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, mcarbon_pos), mcarbon_logic, "mid_carbon_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Middle Carbon pos :" << mcarbon_pos << std::endl;

  //Middle Mylar
  G4Tubs *mmylar_solid = new G4Tubs("middle_mylar_solid", Rin, Rout, 0.5 * mylar_thick, 0., 360 * deg);
  G4LogicalVolume *mmylar_logic = new G4LogicalVolume(mmylar_solid, G4_Mylar, "middle_Mylar_Logic", 0, 0, 0);
  mmylar_logic->SetVisAttributes(mylar_att);
  double mmylar_pos = mcarbon_pos + 0.5 * carbon_thick + 0.5 * mylar_thick;

  mmylar_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, mmylar_pos), mmylar_logic, "middle_mylar_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Middle mylar pos :" << mmylar_pos << std::endl;

  //Middle PCB front Cu
  G4Tubs *mpcbcu_solid = new G4Tubs("middle_pcbcu_solid", Rin, Rout, 0.5 * cu_thick, 0., 360 * deg);
  G4LogicalVolume *mpcbcu_logic = new G4LogicalVolume(mpcbcu_solid, G4_Cu, "middle_PCBCu_Logic", 0, 0, 0);
  mpcbcu_logic->SetVisAttributes(cu_att);
  double mpcbcu_pos = mmylar_pos + 0.5 * mylar_thick + 0.5 * cu_thick;

  mpcbcu_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, mpcbcu_pos), mpcbcu_logic, "middle_pcbcu_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Middle PCB Cu pos :" << mpcbcu_pos << std::endl;

  //Middle PCB
  G4Tubs *mpcb_solid = new G4Tubs("middle_pcb_solid", Rin, Rout, 0.5 * pcb_thick, 0., 360 * deg);
  G4LogicalVolume *mpcb_logic = new G4LogicalVolume(mpcb_solid, G4_FR4, "middle_PCB_Logic", 0, 0, 0);
  mpcb_logic->SetVisAttributes(tof_pcb_att);
  double mpcb_pos = mpcbcu_pos + 0.5 * cu_thick + 0.5 * pcb_thick;

  mpcb_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, mpcb_pos), mpcb_logic, "middle_pcb_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Middle PCB pos :" << mpcb_pos << std::endl;

  //Middle PCB back  Cu
  G4Tubs *mpcbcu2_solid = new G4Tubs("middle_pcbcu2_solid", Rin, Rout, 0.5 * cu_thick, 0., 360 * deg);
  G4LogicalVolume *mpcbcu2_logic = new G4LogicalVolume(mpcbcu2_solid, G4_Cu, "middle_Pcbcu2_Logic", 0, 0, 0);
  mpcbcu2_logic->SetVisAttributes(cu_att);
  double mpcbcu2_pos = mpcb_pos + 0.5 * pcb_thick + 0.5 * cu_thick;

  mpcbcu2_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, mpcbcu2_pos), mpcbcu2_logic, "middle_pcbcu2_phys", Logic, 0, false, OverlapCheck());
  if (Verbosity() > 1) std::cout << " Middle PCB Cu pos :" << mpcbcu2_pos << std::endl;

  G4Tubs *mmylar2_solid = new G4Tubs("middle_mylar2_solid", Rin, Rout, 0.5 * mylar_thick, 0., 360 * deg);
  G4LogicalVolume *mmylar2_logic = new G4LogicalVolume(mmylar2_solid, G4_Mylar, "middle_Mylar2_Logic", 0, 0, 0);
  mmylar2_logic->SetVisAttributes(mylar_att);
  double mmylar2_pos = mpcbcu2_pos + 0.5 * cu_thick + 0.5 * mylar_thick;

  mmylar2_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, mmylar2_pos), mmylar2_logic, "middle_mylar2_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Middle mylar pos :" << mmylar2_pos << std::endl;

  //Middle PCB carbon
  G4Tubs *mcarbon2_solid = new G4Tubs("mid_carbon2_solid", Rin, Rout, 0.5 * carbon_thick, 0., 360 * deg);
  G4LogicalVolume *mcarbon2_logic = new G4LogicalVolume(mcarbon2_solid, G4_Ca, "mid_Carbon2_Logic", 0, 0, 0);
  mcarbon2_logic->SetVisAttributes(carbon_att);
  double mcarbon2_pos = mmylar2_pos + 0.5 * mylar_thick + 0.5 * carbon_thick;
  mcarbon2_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, mcarbon2_pos), mcarbon2_logic, "mid_carbon2_phys", Logic, 0, false, OverlapCheck());
  if (Verbosity() > 1) std::cout << " Middle PCB back Carbon pos :" << mcarbon2_pos << std::endl;

  //Back gas and mrpc layers
  double bglass_begin_pos = mcarbon2_pos + 0.5 * carbon_thick;
  double bglass_mid_pos = 0.;
  double bgas_begin_pos = bglass_begin_pos + glass_thick;
  double bgas_mid_pos = 0.;

  if (Verbosity() > 1)
  {
    std::cout << " Back glass begin pos : " << bglass_begin_pos << std::endl;
    std::cout << " Back gas begin pos : " << bgas_begin_pos << std::endl;
  }

  G4Tubs *bglass_solid[n_bgas_layer + 1];
  G4LogicalVolume *bglass_logic[n_bgas_layer + 1];
  G4Tubs *bgas_solid[n_bgas_layer];
  G4LogicalVolume *bgas_logic[n_bgas_layer];

  for (int l = 0; l < n_bgas_layer + 1; l++)
  {
    /*
      bglass_solid[l] = new G4Tubs(sprintf("back_glass_solid_%d", l), Rin, Rout, 0.5*glass_thick, 0, 360.*deg);
      bglass_logic[l] = new G4LogicalVolume(bglass_solid[l], plateglass, sprintf("back_glass_Logic_%d", l),0 , 0,0);
    
      bglass_solid[l] = new G4Tubs("back_glass_solid", Rin, Rout, 0.5*glass_thick, 0, 360.*deg);
      bglass_logic[l] = new G4LogicalVolume(bglass_solid[l], plateglass, "back_glass_Logic",0 , 0,0);
    */
    std::string bgl_solid_name = "bglass_solid_" + std::to_string(l);
    std::string bgl_logic_name = "bglass_logic_" + std::to_string(l);
    std::string bgl_phys_name = "bglass_phys_" + std::to_string(l);

    bglass_solid[l] = new G4Tubs(bgl_solid_name, Rin, Rout, 0.5 * glass_thick, 0, 360. * deg);
    bglass_logic[l] = new G4LogicalVolume(bglass_solid[l], plateglass, bgl_logic_name, 0, 0, 0);
    bglass_logic[l]->SetVisAttributes(glass_att);
    bglass_mid_pos = bglass_begin_pos + (l + 0.5) * glass_thick + l * gas_gap;
    bglass_phys[l] = new G4PVPlacement(0, G4ThreeVector(0, 0, bglass_mid_pos), bglass_logic[l], bgl_phys_name, Logic, 0, false, OverlapCheck());

    if (l < n_bgas_layer)
    {
      /*
 	bgas_solid[l] = new G4Tubs(sprintf("front_gas_solid_%d", l), Rin, Rout, 0.5*gas_gap, 0, 360.*deg);
 	bgas_logic[l] = new G4LogicalVolume(bgas_solid[l], G4_gas, sprintf("back_gas_Logic_%d", l),0 , 0,0);
 		      
        bgas_solid[l] = new G4Tubs("front_gas_solid", Rin, Rout, 0.5*gas_gap, 0, 360.*deg);
        bgas_logic[l] = new G4LogicalVolume(bgas_solid[l], G4_gas, "back_gas_Logic",0 , 0,0);
      */
      std::string bgas_solid_name = "bgas_solid_" + std::to_string(l);
      std::string bgas_logic_name = "bgas_logic_" + std::to_string(l);
      std::string bgas_phys_name = "bgas_phys_" + std::to_string(l);

      bgas_solid[l] = new G4Tubs(bgas_solid_name, Rin, Rout, 0.5 * gas_gap, 0, 360. * deg);
      bgas_logic[l] = new G4LogicalVolume(bgas_solid[l], G4_gas, bgas_logic_name, 0, 0, 0);
      bgas_logic[l]->SetVisAttributes(gas_att);
      bgas_mid_pos = bgas_begin_pos + (l + 0.5) * gas_gap + l * glass_thick;
      // bgas_phys[l] =  new G4PVPlacement(0, G4ThreeVector(0, 0, bgas_mid_pos), bgas_logic[l], sprintf("Gasback_layer_phys_%d",l), Logic, 0, false, OverlapCheck());
      bgas_phys[l] = new G4PVPlacement(0, G4ThreeVector(0, 0, bgas_mid_pos), bgas_logic[l], bgas_phys_name, Logic, 0, false, OverlapCheck());
    }
    if (Verbosity() > 1)
    {
      std::cout << " Back glass layer : " << l << " glass mid : " << bglass_mid_pos << std::endl;
      std::cout << "Back gas layer : " << l << " gas mid : " << bgas_mid_pos << std::endl;
    }
  }

  //Back carbon layer
  G4Tubs *bcarbon_solid = new G4Tubs("mid_carbon2_solid", Rin, Rout, 0.5 * carbon_thick, 0., 360 * deg);
  G4LogicalVolume *bcarbon_logic = new G4LogicalVolume(bcarbon_solid, G4_Ca, "bach_Carbon_Logic", 0, 0, 0);
  bcarbon_logic->SetVisAttributes(carbon_att);
  double bcarbon_pos = bglass_mid_pos + 0.5 * glass_thick + 0.5 * carbon_thick;
  bcarbon_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, bcarbon_pos), bcarbon_logic, "back_carbon_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " back carbon mid :" << bcarbon_pos << std::endl;

  //Back mylar layer
  G4Tubs *bmylar_solid = new G4Tubs("back_mylar_solid", Rin, Rout, 0.5 * mylar_thick, 0., 360 * deg);
  G4LogicalVolume *bmylar_logic = new G4LogicalVolume(bmylar_solid, G4_Mylar, "back_Mylar_Logic", 0, 0, 0);
  bmylar_logic->SetVisAttributes(mylar_att);
  double bmylar_pos = bcarbon_pos + 0.5 * carbon_thick + 0.5 * mylar_thick;

  bmylar_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, bmylar_pos), bmylar_logic, "back_mylar_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Back mylar pos :" << bmylar_pos << std::endl;

  //Back PCB Cu layer
  G4Tubs *bpcbcu_solid = new G4Tubs("back_pcbcu_solid", Rin, Rout, 0.5 * cu_thick, 0., 360 * deg);
  G4LogicalVolume *bpcbcu_logic = new G4LogicalVolume(bpcbcu_solid, G4_Cu, "back_PCBCu_Logic", 0, 0, 0);
  bpcbcu_logic->SetVisAttributes(cu_att);
  double bpcbcu_pos = bmylar_pos + 0.5 * mylar_thick + 0.5 * cu_thick;

  bpcbcu_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, bpcbcu_pos), bpcbcu_logic, "middle_pcbcu_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Middle PCB back Cu pos :" << bpcbcu_pos << std::endl;

  //Back PCB
  G4Tubs *bpcb_solid = new G4Tubs("back_pcb_solid", Rin, Rout, 0.5 * pcb_thick, 0., 360 * deg);
  G4LogicalVolume *bpcb_logic = new G4LogicalVolume(bpcb_solid, G4_FR4, "back_PCB_Logic", 0, 0, 0);
  bpcb_logic->SetVisAttributes(tof_pcb_att);
  double bpcb_pos = bpcbcu_pos + 0.5 * cu_thick + 0.5 * pcb_thick;

  bpcb_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, bpcb_pos), bpcb_logic, "back_pcb_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Back PCB pos :" << bpcb_pos << std::endl;

  //Back Honeycomb
  G4Tubs *bhc_solid = new G4Tubs("front_honeycomb_solid", Rin, Rout, 0.5 * honeycomb_thick, 0., 360 * deg);
  G4LogicalVolume *bhc_logic = new G4LogicalVolume(bhc_solid, G4_Nomex, "front_honeycomb_Logic", 0, 0, 0);
  bhc_logic->SetVisAttributes(tof_hc_att);
  double bhc_pos = bpcb_pos + 0.5 * pcb_thick + 0.5 * honeycomb_thick;

  bhc_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, bhc_pos), bhc_logic, "front_hc_phys", Logic, 0, false, OverlapCheck());

  if (Verbosity() > 1) std::cout << " Front HC mid pos :" << bhc_pos << std::endl;
}
