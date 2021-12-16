#include "PHG4BarrelEcalDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

PHG4BarrelEcalDisplayAction::PHG4BarrelEcalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4BarrelEcalDisplayAction::~PHG4BarrelEcalDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4BarrelEcalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    if (it.second == "BCalCylinder")
    {
      visatt->SetColour(G4Colour::White());
      visatt->SetForceWireframe(true);
      visatt->SetVisibility(false);
    }
    else if (it.second == "Block1")
    {
      visatt->SetColour(G4Colour(1., 0.5, 0.));
    }
    else if (it.second == "Block2")
    {
      visatt->SetColour(G4Colour::Blue());
    }
    else if (it.second == "Glass")
    {
      visatt->SetColour(G4Colour::White());
      //visatt->SetForceWireframe(true);
    }
    else if (it.second == "Si")
    {
      visatt->SetColour(G4Colour::Green());
    }
    else if (it.second == "Kapton")
    {
      visatt->SetColour(G4Colour::Brown());
    }
    else if (it.second == "SIO2")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "Carbon")
    {
      visatt->SetColour(G4Colour::Black());
    }
    else
    {
      std::cout << "unknown logical volume " << it.second << std::endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
