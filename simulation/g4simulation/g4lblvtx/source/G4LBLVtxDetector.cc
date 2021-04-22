#include "G4LBLVtxDetector.h"

#include "G4LBLVtxDisplayAction.h"

#include <g4main/PHG4Detector.h>  // for PHG4Detector
#include <g4main/PHG4Subsystem.h>
#include "g4main/PHG4DisplayAction.h"  // for PHG4DisplayAction

#include <phparameter/PHParameters.h>

#include <TSystem.h>

#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <iostream>  // for operator<<, basic_ostream
#include <memory>

using namespace std;

G4LBLVtxDetector::G4LBLVtxDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam, PHParameters* parameters)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<G4LBLVtxDisplayAction*>(subsys->GetDisplayAction()))
  , m_GDMPath(parameters->get_string_param("GDMPath"))
  , m_TopVolName(parameters->get_string_param("TopVolName"))
  , m_placeX(parameters->get_double_param("place_x") * cm)
  , m_placeY(parameters->get_double_param("place_y") * cm)
  , m_placeZ(parameters->get_double_param("place_z") * cm)
  , m_rotationX(parameters->get_double_param("rot_x") * degree)
  , m_rotationY(parameters->get_double_param("rot_y") * degree)
  , m_rotationZ(parameters->get_double_param("rot_z") * degree)
  , m_Active(parameters->get_int_param("active"))
  , m_AbsorberActive(parameters->get_int_param("absorberactive"))
{
  m_ActiveVolName.insert("TopFilament");
  m_ActiveVolName.insert("CarbonFleeceBottom");
}

G4LBLVtxDetector::~G4LBLVtxDetector()
{
}

void G4LBLVtxDetector::Print(const std::string& what) const
{
  cout << "G4LBLVtxDetector::" << GetName() << " - import " << m_TopVolName << " from " << m_GDMPath << " with shift "
       << m_placeX << ","
       << m_placeY << ","
       << m_placeZ << "cm and rotation "
       << m_rotationX << ","
       << m_rotationY << ","
       << m_rotationZ << "rad" << endl;
}

void G4LBLVtxDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << " G4LBLVtxDetector::Construct:";
    Print();
  }

  //===================================
  // Import the stave physical volume here
  //===================================

  // import the staves from the geometry file
  unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.SetOverlapCheck(OverlapCheck());
  gdmlParser.Read(m_GDMPath, OverlapCheck());

  G4LogicalVolume* vol = reader->GetVolume(m_TopVolName);

  if (!vol)
  {
    cout << "G4LBLVtxDetector::Construct - Fatal Error - failed to find G4LogicalVolume " << m_TopVolName << " - Print: ";
    Print();
    gSystem->Exit(1);
  }
  PHG4Subsystem* mysys = GetMySubsystem();
  mysys->SetLogicalVolume(vol);

  G4RotationMatrix* rotm = new G4RotationMatrix();
  rotm->rotateX(m_rotationX);
  rotm->rotateY(m_rotationY);
  rotm->rotateZ(m_rotationZ);
  G4ThreeVector placeVec(m_placeX, m_placeY, m_placeZ);

  G4VPhysicalVolume* phys = new G4PVPlacement(rotm, placeVec,
                                              vol,
                                              G4String(GetName().c_str()),
                                              logicWorld, false, 0, OverlapCheck());
  SetActiveVolumes(phys);
}

void G4LBLVtxDetector::SetActiveVolumes(G4VPhysicalVolume* physvol)
{
  G4LogicalVolume* logvol = physvol->GetLogicalVolume();
  if (logvol->GetNoDaughters() == 0)  // add only if no other volumes inside
  {
    m_DisplayAction->AddLogVolume(logvol);
    string test(physvol->GetName());
    int added_to_active = 0;
    for (set<std::string>::const_iterator iter = m_ActiveVolName.begin();
         iter != m_ActiveVolName.end(); ++iter)
    {
      //      cout << "checking " << test << " for " << *iter << endl;
      if (test.find(*iter) != string::npos)
      {
        //	cout << "Adding " << physvol->GetName() << endl;
        m_ActivePhysVolumeMap.insert(physvol);
        added_to_active = 1;
      }
    }
    if (!added_to_active)
    {
      m_PassivePhysVolumeMap.insert(physvol);
    }
  }
  else
  {
    for (unsigned int i = 0; i < (unsigned int) logvol->GetNoDaughters(); ++i)
    {
      G4VPhysicalVolume* physvol = logvol->GetDaughter(i);
      // here we decide which volumes are active
      SetActiveVolumes(physvol);
    }
  }
  return;
}

int G4LBLVtxDetector::IsInDetector(G4VPhysicalVolume* physvol) const
{
  if (m_Active)
  {
    if (m_ActivePhysVolumeMap.find(physvol) != m_ActivePhysVolumeMap.end())
    {
      return 1;
    }
  }
  if (m_AbsorberActive)
  {
    if (m_PassivePhysVolumeMap.find(physvol) != m_PassivePhysVolumeMap.end())
    {
      return -1;
    }
  }
  return 0;
}
