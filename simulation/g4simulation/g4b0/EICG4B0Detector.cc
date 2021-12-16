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

#include "EICG4B0Detector.h"
//#include "EICG4B0SteppingAction.h"
//#include "EICG4B0Subsystem.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
//#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Color.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <TSystem.h>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <utility>
 
class G4VSolid;
class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
EICG4B0Detector::EICG4B0Detector(PHG4Subsystem *subsys,
                                 PHCompositeNode *Node,
                                 PHParameters *parameters,
                                 const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_Layer(lyr)
{
}

//_______________________________________________________________
int EICG4B0Detector::IsInDetector(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end())
  {
    return 1;
  }
  return 0;
}

int EICG4B0Detector::GetDetId(G4VPhysicalVolume *volume) const
{
  if (IsInDetector(volume))
  {
    return 1;
  }
  return -1;
}

//_______________________________________________________________
void EICG4B0Detector::ConstructMe(G4LogicalVolume *logicWorld)
{
  //begin implement your own here://
  // Do not forget to multiply the parameters with their respective CLHEP/G4 unit !
  if (Verbosity() > 0)
  {
    std::cout << "EICG4B0Detector: Begin Construction" << std::endl;
  }

  cout << " !!! length = " << m_Params->get_double_param("length");
  if (m_Params->get_double_param("spanningAngle") >= 360)
  {
    cout << " !!! No PACKMAN" << endl;
    return;
  }
  //Print("ALL");


  G4VSolid *solid0 = new G4Tubs("EICG4B0Solid0",
                                0.,
                                m_Params->get_double_param("outer_radius") * cm,
                                m_Params->get_double_param("length") / 2. * cm,
                                m_Params->get_double_param("startAngle") * degree,
                                m_Params->get_double_param("spanningAngle") * degree);
  G4VSolid *solidPipeHole = new G4Tubs("EICG4B0IonPipeSolid",
                                       0.,
                                       m_Params->get_double_param("pipe_hole_r") * cm,
                                       m_Params->get_double_param("length") * cm,
                                       0., 360. * degree);
  G4VSolid *solidCableHole = new G4Tubs("EICG4B0CableSolid",
                                       0.,
                                       m_Params->get_double_param("cable_hole") * cm,
                                       m_Params->get_double_param("length") * cm,
                                       0., 360. * degree);
  G4VSolid *solidPipeHole1 = new G4Box("EICG4B0PipeSolid1",
                                m_Params->get_double_param("pipe_hole") / 2. * cm,
                                m_Params->get_double_param("pipe_hole_r") * cm,
                                m_Params->get_double_param("length") * cm);
  G4VSolid *solid1 = new G4Tubs("EICG4B0Solid1",
                                0.,
                                (m_Params->get_double_param("outer_radius") - m_Params->get_double_param("d_radius")) * cm,
                                m_Params->get_double_param("length") / 2. * cm,
                                (m_Params->get_double_param("startAngle") + m_Params->get_double_param("spanningAngle")) * degree,
                                (360 - m_Params->get_double_param("spanningAngle")) * degree);
  G4UnionSolid *solid10 = new G4UnionSolid("EICG4B0Solid10", solid0, solid1);
//here we subtract two cilinders and one box from the detector volume to get an oval-shaped hole
  G4SubtractionSolid *solids = new G4SubtractionSolid("EICG4B0Solid", solid10, solidPipeHole, 0, G4ThreeVector((m_Params->get_double_param("pipe_x")+m_Params->get_double_param("pipe_hole")/2) * cm, m_Params->get_double_param("pipe_y") * cm, m_Params->get_double_param("pipe_z") * cm));
  G4SubtractionSolid *solids1 = new G4SubtractionSolid("EICG4B0Solid", solids, solidPipeHole, 0, G4ThreeVector((m_Params->get_double_param("pipe_x")-m_Params->get_double_param("pipe_hole")/2) * cm, m_Params->get_double_param("pipe_y") * cm, m_Params->get_double_param("pipe_z") * cm));
  G4SubtractionSolid *solids2 = new G4SubtractionSolid("EICG4B0Solid", solids1, solidCableHole, 0, G4ThreeVector(m_Params->get_double_param("cable_x") * cm, m_Params->get_double_param("cable_y") * cm, m_Params->get_double_param("cable_z") * cm));
  G4SubtractionSolid *solidB0 = new G4SubtractionSolid("EICG4B0Solid", solids2, solidPipeHole1, 0, G4ThreeVector(m_Params->get_double_param("pipe_x") * cm, m_Params->get_double_param("pipe_y") * cm, m_Params->get_double_param("pipe_z") * cm));
  G4RotationMatrix *rotm = new G4RotationMatrix();
  G4LogicalVolume *logical = new G4LogicalVolume(solidB0,
                                                 G4Material::GetMaterial(m_Params->get_string_param("material")),
                                                 "EICG4B0Logical");

  G4VisAttributes *vis = new G4VisAttributes(G4Color(0.8, 0.4, 0.2, 1.0));
  if (m_Params->get_string_param("material") == "G4_PbWO4") vis->SetColor(0.8, 0.4, 0.2, 1.0);
  if (m_Params->get_string_param("material") == "G4_Cu") vis->SetColor(1., 0., 1., .5);
  if (m_Params->get_string_param("material") == "G4_Si") vis->SetColor(1., 1., 0., .8);
  vis->SetForceSolid(true);
  logical->SetVisAttributes(vis);
  rotm->rotateY(m_Params->get_double_param("rot_y") * deg);
  //  rotm->rotateY(0. * deg);

  G4VPhysicalVolume *phy = new G4PVPlacement(
      rotm,
      G4ThreeVector(m_Params->get_double_param("place_x") * cm,
                    m_Params->get_double_param("place_y") * cm,
                    m_Params->get_double_param("place_z") * cm),
      //      G4ThreeVector(0. * cm,
      //                    0. * cm,
      //                    0. * cm),
      logical, "EICG4B0", logicWorld, 0, false, OverlapCheck());
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

//_______________________________________________________________
void EICG4B0Detector::Print(const std::string &what) const
{
  std::cout << "EICG4B0 Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}

PHParameters *EICG4B0Detector::getParams()
{
  return m_Params;
}
