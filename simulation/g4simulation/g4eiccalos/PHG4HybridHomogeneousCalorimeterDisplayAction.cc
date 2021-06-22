#include "PHG4HybridHomogeneousCalorimeterDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

using namespace std;

PHG4HybridHomogeneousCalorimeterDisplayAction::PHG4HybridHomogeneousCalorimeterDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4HybridHomogeneousCalorimeterDisplayAction::~PHG4HybridHomogeneousCalorimeterDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4HybridHomogeneousCalorimeterDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    if (it.second == "CarbonShell")
    {
      visatt->SetColour(G4Colour::Black());
    }
    else if (it.second == "Crystal")
    {
      visatt->SetColour(G4Colour::Cyan());
    }
    else if (it.second == "Envelope")
    {
      visatt->SetVisibility(false);
    }
    else if (it.second == "Invisible")
    {
      visatt->SetVisibility(false);
    }
    else if (it.second == "TwoByTwo")
    {
      visatt->SetColour(G4Colour::Gray());
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
