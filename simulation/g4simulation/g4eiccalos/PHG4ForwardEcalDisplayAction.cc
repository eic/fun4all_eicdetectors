#include "PHG4ForwardEcalDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

PHG4ForwardEcalDisplayAction::PHG4ForwardEcalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4ForwardEcalDisplayAction::PHG4ForwardEcalDisplayAction(const std::string &name, bool detailed)
  : PHG4DisplayAction(name)
{
  showdetails = detailed;
  if (!detailed) std::cout << "PHG4ForwardEcalDisplayAction::disabled detailed view of towers" << std::endl;
}


PHG4ForwardEcalDisplayAction::~PHG4ForwardEcalDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4ForwardEcalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    if (it.second == "Absorber")
    {
      visatt->SetColour(G4Colour::Gray());
      if (showdetails)
        visatt->SetVisibility(true);
      else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "Coating")
    {
      visatt->SetColour(G4Colour::Black());
      if (showdetails)
        visatt->SetVisibility(true);
      else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "Clamp")
    {
      visatt->SetColour(4 * 21. / 255, 4 * 27. / 255, 4 * 31. / 255);
      // visatt->SetForceWireframe(true);
      if (showdetails)
        visatt->SetVisibility(true);
      else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "Envelope")
    {
      visatt->SetColour(G4Colour::Magenta());
      visatt->SetVisibility(false);
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "Fiber")
    {
      // visatt->SetColour(G4Colour::Cyan());
      visatt->SetColour(152./ 255,251./ 255,152./ 255, 0.4);
      if (showdetails)
        visatt->SetVisibility(true);
      else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "Scintillator")
    {
      // visatt->SetColour(G4Colour::White());
      // visatt->SetColour(G4Colour::Cyan());
      visatt->SetColour(127./ 255,255./ 255,212./ 255, 0.2);
      if (showdetails)
        visatt->SetVisibility(true);
      else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "SingleTower")
    {
      // visatt->SetColour(G4Colour::Cyan());
      visatt->SetColour(4 * 21. / 255, 4 * 27. / 255, 4 * 31. / 255);
      // visatt->SetForceWireframe(true);
      if (showdetails)
        visatt->SetVisibility(false);
      else 
        visatt->SetVisibility(true);
    }
    else if (it.second == "miniblock")
    {
      visatt->SetVisibility(false);
      visatt->SetColour(G4Colour::Red());
      visatt->SetForceSolid(false);
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
