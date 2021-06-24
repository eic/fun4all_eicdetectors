#include "PHG4TTLDetector.h"
#include "PHG4TTLDisplayAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>      // for PHG4Subsystem

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

#include <map>  // for _Rb_tree_iterator, _Rb_tree_co...
#include <utility>  // for pair
#include <algorithm>  // for max
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <climits>

class G4VPhysicalVolume;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4TTLDetector::PHG4TTLDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters* parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , name_base(dnam)
  , m_DisplayAction(dynamic_cast<PHG4TTLDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , overlapcheck_sector(false)
{
}

//_______________________________________________________________
//_______________________________________________________________
bool PHG4TTLDetector::IsInSectorActive(G4VPhysicalVolume *physvol)
{
  for (map_phy_vol_t::const_iterator it = map_active_phy_vol.begin();
       it != map_active_phy_vol.end(); ++it)
  {
    if (physvol == (*it).second)
    {
      return true;
    }
  }

  return false;
}

//_______________________________________________________________
void PHG4TTLDetector::ConstructMe(G4LogicalVolume *logicWorld)
{

  const int nLayers = 7;
  string strLayerName[nLayers] = {
    "SiliconSensor", "Metalconnection", "HDI", "Cooling", "Support", "Support_Gap", "Support2"
  };
  G4Material* materialLayer[nLayers] = {
    G4Material::GetMaterial("G4_Si"), G4Material::GetMaterial("G4_Al"), G4Material::GetMaterial("G4_KAPTON"), G4Material::GetMaterial("G4_WATER"), G4Material::GetMaterial("G4_GRAPHITE"), G4Material::GetMaterial("G4_AIR"), G4Material::GetMaterial("G4_GRAPHITE")
  };
  G4double thicknessLayer[nLayers] = {
    m_Params->get_double_param("tSilicon"), 100 * um, 20 * um, 100 * um, 50 * um, 1 * cm, 50 * um
  };
  bool layerActive[nLayers] = {
    true, false, false, false, false, false, false
  };

  G4double thicknessDet = 0;
  for(int ilay=0; ilay<nLayers; ilay++){
    thicknessDet += thicknessLayer[ilay]/10;
  }

  //Create the envelope = 'world volume' for the calorimeter
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4double rMin = m_Params->get_double_param("rMin");
  G4double rMax = m_Params->get_double_param("rMax");
  G4double xoffset = m_Params->get_double_param("offset_x");

  G4VSolid* beampipe_cutout = new G4Cons("ttl_beampipe_cutout",
                                      0, rMin,
                                      0, rMin,
                                      thicknessDet / 2.0,
                                      0, 2 * M_PI);
  G4VSolid* ttl_envelope_solid = new G4Cons("ttl_beampipe_cutout",
                                      0, rMax,
                                      0, rMax,
                                      thicknessDet / 2.0,
                                      0, 2 * M_PI);
  ttl_envelope_solid = new G4SubtractionSolid(G4String("ttl_envelope_solid"), ttl_envelope_solid, beampipe_cutout
                                                            , 0 ,G4ThreeVector( xoffset , 0 ,0.));

  const G4Transform3D transform_Det_to_Hall =
      G4RotateX3D(-m_Params->get_double_param("polar_angle")) * G4TranslateZ3D(
                                                        m_Params->get_double_param("place_z") + thicknessDet / 2);

  G4LogicalVolume *DetectorLog_Det = new G4LogicalVolume(ttl_envelope_solid,  //
                                                         Air, name_base + "_Log");
  RegisterLogicalVolume(DetectorLog_Det);

  RegisterPhysicalVolume(
        new G4PVPlacement(
            G4RotateZ3D(2 * pi) * transform_Det_to_Hall, DetectorLog_Det,
            name_base + "_Physical", logicWorld, false, 0, overlapcheck_sector));

  double z_start = -thicknessDet / 2;
  for(int ilay=0; ilay<nLayers; ilay++){
    const string layer_name = name_base + "_" + strLayerName[ilay];
    const string layer_name_Solid = layer_name + "_Sol";

    G4VSolid *Sol_Raw = new G4Tubs(layer_name_Solid + "_Raw",     //const G4String& pName,
                                  0,                 //      G4double pRMin,
                                  rMax,  //      G4double pRMax,
                                  thicknessLayer[ilay] / 10 / 2,     //      G4double pDz,
                                  0,                 //      G4double pSPhi,
                                  2 * pi             //      G4double pDPhi
    );

    G4VSolid *LayerSol_Det = new G4SubtractionSolid(layer_name_Solid.c_str(), Sol_Raw, beampipe_cutout
                                                            , 0 ,G4ThreeVector( xoffset , 0 ,0.));
    G4LogicalVolume *LayerLog_Det = new G4LogicalVolume(LayerSol_Det,  //
                                                        materialLayer[ilay], layer_name + "_Log");
    RegisterLogicalVolume(LayerLog_Det);
    RegisterPhysicalVolume(
        new G4PVPlacement(0, G4ThreeVector(0, 0, z_start + thicknessLayer[ilay] / 10 / 2), LayerLog_Det,
                          layer_name + "_Physical", DetectorLog_Det, false, 0, overlapcheck_sector),
        layerActive[ilay]);
    z_start += thicknessLayer[ilay]/10;
  }

  m_DisplayAction->AddVolume(DetectorLog_Det, "DetectorBox");

  for (map_log_vol_t::iterator it = map_log_vol.begin(); it != map_log_vol.end(); ++it)
  {
    if ((*it).first != G4String(name_base + "_Log"))
    {
      // sub layers
      m_DisplayAction->AddVolume((*it).second, "TTLDetector");
    }
  }
}



G4LogicalVolume *
PHG4TTLDetector::RegisterLogicalVolume(G4LogicalVolume *v)
{
  if (!v)
  {
    std::cout
        << "PHG4TTLDetector::RegisterVolume - Error - invalid volume!"
        << std::endl;
    return v;
  }
  if (map_log_vol.find(v->GetName()) != map_log_vol.end())
  {
    std::cout << "PHG4TTLDetector::RegisterVolume - Warning - replacing " << v->GetName() << std::endl;
  }

  map_log_vol[v->GetName()] = v;

  return v;
}

G4PVPlacement *
PHG4TTLDetector::RegisterPhysicalVolume(G4PVPlacement *v,
                                              const bool active)
{
  if (!v)
  {
    std::cout
        << "PHG4TTLDetector::RegisterPhysicalVolume - Error - invalid volume!"
        << std::endl;
    return v;
  }

  phy_vol_idx_t id(v->GetName(), v->GetCopyNo());

  if (map_phy_vol.find(id) != map_phy_vol.end())
  {
    std::cout
        << "PHG4TTLDetector::RegisterPhysicalVolume - Warning - replacing "
        << v->GetName() << "[" << v->GetCopyNo() << "]" << std::endl;
  }

  map_phy_vol[id] = v;

  if (active)
    map_active_phy_vol[id] = v;

  return v;
}
