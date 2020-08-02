#include "BeastMagnetDetector.h"

#include "BeastMagnetDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4Subsystem.h>

#include <TSystem.h>
#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4LogicalVolumeStore.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SolidStore.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>
#include <memory>

class G4VSolid;
class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
BeastMagnetDetector::BeastMagnetDetector(PHG4Subsystem *subsys,
                                         PHCompositeNode *Node,
                                         PHParameters *parameters,
                                         const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<BeastMagnetDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_GDMPath(parameters->get_string_param("GDMPath"))
  , m_TopVolName(parameters->get_string_param("TopVolName"))
  , m_AbsorberActive(parameters->get_int_param("absorberactive"))
{
}

//_______________________________________________________________
int BeastMagnetDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end())
  {
    return 1;
  }
  return 0;
}

//_______________________________________________________________
void BeastMagnetDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.SetOverlapCheck(OverlapCheck());
  gdmlParser.Read(m_GDMPath, false);

  // alright the reader just puts everything into G4 Stores - endless fun
  // print the show out, first solids:
  // for (auto i=G4SolidStore::GetInstance()->begin(); i!=G4SolidStore::GetInstance()->end(); ++i)
  //   cout << "solid vol name: " << (*i)->GetName() << endl;
  // for (auto i=G4LogicalVolumeStore::GetInstance()->begin(); i!=G4LogicalVolumeStore::GetInstance()->end(); i++)
  //   cout << "logvol name " << (*i)->GetName() << endl;

  G4AssemblyVolume *avol = reader->GetAssembly(m_TopVolName);
  if (!avol)
  {
    cout << "not found" << endl;
  }

  G4RotationMatrix *Rot = new G4RotationMatrix();
  G4ThreeVector g4vec;
  avol->MakeImprint(logicWorld, g4vec, Rot, 0, OverlapCheck());
  vector<G4VPhysicalVolume *>::iterator it = avol->GetVolumesIterator();
  for (unsigned int i = 0; i < avol->TotalImprintedVolumes(); i++)
  {
    InsertVolumes(*it);
    ++it;
  }

  return;
}

void BeastMagnetDetector::InsertVolumes(G4VPhysicalVolume *physvol)
{
  G4LogicalVolume *logvol = physvol->GetLogicalVolume();
  m_DisplayAction->AddLogicalVolume(logvol);
  m_PhysicalVolumesSet.insert(physvol);
  // G4 10.06 returns unsigned int for GetNoDaughters()
  // lower version return int, need to cast to avoid compiler error
  for (int i = 0; i < (int) logvol->GetNoDaughters(); ++i)
  {
    G4VPhysicalVolume *physvol = logvol->GetDaughter(i);
    // here we decide which volumes are active
    InsertVolumes(physvol);
  }
  return;
}

//_______________________________________________________________
void BeastMagnetDetector::Print(const std::string &what) const
{
  cout << "BeastMagnet Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Version 0.1" << endl;
    cout << "Parameters:" << endl;
    m_Params->Print();
  }
  return;
}
