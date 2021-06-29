#include "PHG4ForwardHcalDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

PHG4ForwardHcalDisplayAction::PHG4ForwardHcalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4ForwardHcalDisplayAction::PHG4ForwardHcalDisplayAction(const std::string &name, bool detailed)
  : PHG4DisplayAction(name)
{
  showdetails = detailed;
  if (!detailed) std::cout << "PHG4ForwardHcalDisplayAction::disabled detailed view of towers" << std::endl;
}

PHG4ForwardHcalDisplayAction::~PHG4ForwardHcalDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4ForwardHcalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    else if (it.second == "FHcalEnvelope")
    {
      visatt->SetVisibility(false);
    }
    else if (it.second == "Scintillator")
    {
      visatt->SetColour(G4Colour::White());
      if (showdetails)
        visatt->SetVisibility(true);
      else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "WLSplate")
    {
      visatt->SetColour(G4Colour::Yellow());
      if (showdetails)
        visatt->SetVisibility(true);
      else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "SupportPlate")
    {
      visatt->SetColour(G4Colour::Gray());
      if (showdetails)
        visatt->SetVisibility(true);
      else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "SingleTower")
    {
      visatt->SetColour(G4Colour::Gray());
      if (showdetails)
        visatt->SetVisibility(false);
      else 
        visatt->SetVisibility(true);
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
