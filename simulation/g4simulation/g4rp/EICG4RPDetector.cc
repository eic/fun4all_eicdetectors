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
// Parameter names and defaults are set in EICG4B0Subsystem::SetDefaultParameters()
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

#include "EICG4RPDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>

class G4VSolid;
class PHCompositeNode;

//____________________________________________________________________________..
EICG4RPDetector::EICG4RPDetector(PHG4Subsystem *subsys,
                                 PHCompositeNode *Node,
                                 PHParameters *parameters,
                                 const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_Layer(lyr)
{
}

//_______________________________________________________________
int EICG4RPDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  std::set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end())
  {
    return 1;
  }
  return 0;
}

int EICG4RPDetector::GetDetId(G4VPhysicalVolume *volume) const
{
  if (IsInDetector(volume))
  {
    return 1;
  }
  return -1;
}

//_______________________________________________________________
void EICG4RPDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  //begin implement your own here://
  // Do not forget to multiply the parameters with their respective CLHEP/G4 unit !
  if (Verbosity() > 1) std::cout << " !!! length = " << m_Params->get_double_param("length");
  G4VSolid *solid0 = new G4Box("EICG4RPSolid0",
                               m_Params->get_double_param("rp_x") / 2. * cm,
                               m_Params->get_double_param("rp_y") / 2. * cm,
                               m_Params->get_double_param("length") / 2. * cm);
  G4VSolid *solidPipeHole = new G4Box("EICG4RPIonPipeSolid",
                                      m_Params->get_double_param("hole_x") / 2. * cm,
                                      m_Params->get_double_param("hole_y") / 2. * cm,
                                      m_Params->get_double_param("length") * cm);
  G4SubtractionSolid *solidRP = new G4SubtractionSolid("EICG4RPSolid", solid0, solidPipeHole, 0, G4ThreeVector(m_Params->get_double_param("pipe_x") * cm, m_Params->get_double_param("pipe_y") * cm, m_Params->get_double_param("pipe_z") * cm));
  G4LogicalVolume *logical = new G4LogicalVolume(solidRP,
                                                 G4Material::GetMaterial(m_Params->get_string_param("material")),
                                                 "EICG4RPLogical");

  G4VisAttributes *vis = new G4VisAttributes(G4Color(0.8, 0.4, 0.2, 1.0));
  if (m_Params->get_string_param("material") == "G4_Cu") vis->SetColor(1., 0., 1., .5);
  if (m_Params->get_string_param("material") == "G4_Si") vis->SetColor(1., 1., 0., .8);
  vis->SetForceSolid(true);
  logical->SetVisAttributes(vis);
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateY(m_Params->get_double_param("rot_y") * deg);

  G4VPhysicalVolume *phy = new G4PVPlacement(
      rotm,
      G4ThreeVector(m_Params->get_double_param("place_x") * cm,
                    m_Params->get_double_param("place_y") * cm,
                    m_Params->get_double_param("place_z") * cm),
      //      G4ThreeVector(0. * cm,
      //                    0. * cm,
      //                    0. * cm),
      logical, "EICG4RP", logicWorld, 0, false, OverlapCheck());
  // add it to the list of placed volumes so the IsInDetector method
  // picks them up
  m_PhysicalVolumesSet.insert(phy);
  // hard code detector id to detid
  m_PhysicalVolumesDet.insert({phy, m_Params->get_double_param("detid") + 1});
  //  m_LogicalVolumesSet.insert(logical);
  //end implement your own here://
  return;
}

//_______________________________________________________________
void EICG4RPDetector::Print(const std::string &what) const
{
  std::cout << "EICG4RP Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}

PHParameters *EICG4RPDetector::getParams()
{
  return m_Params;
}
