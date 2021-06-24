#include "PHG4TTLDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Utils.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>                     // for pair

using namespace std;

PHG4TTLDisplayAction::PHG4TTLDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4TTLDisplayAction::~PHG4TTLDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4TTLDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing
  for (auto it : m_LogicalVolumeMap)
  {
    G4LogicalVolume *logvol = it.first;
    if (logvol->GetVisAttributes())
    {
      continue;
    }
    G4VisAttributes *visatt = new G4VisAttributes();
    visatt->SetVisibility(true);
    visatt->SetForceSolid(true);
    m_VisAttVec.push_back(visatt);  // for later deletion
    if (it.second == "TTLDetector")
    {
      PHG4Utils::SetColour(visatt, it.first->GetMaterial()->GetName());
    }
    else if (it.second == "DetectorBox")
    {
      visatt->SetVisibility(false);
      // visatt->SetColour(G4Colour::White());
      // visatt->SetForceWireframe(true);
      // visatt->SetForceLineSegmentsPerCircle(50);
    }
    else
    {
      cout << "unknown logical volume " << it.second << endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
