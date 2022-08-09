#include "PHG4FSTDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>                     // for pair

using namespace std;

PHG4FSTDisplayAction::PHG4FSTDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4FSTDisplayAction::~PHG4FSTDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4FSTDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    // visatt->SetVisibility(false);
    visatt->SetForceSolid(true);
    m_VisAttVec.push_back(visatt);  // for later deletion
 
  //  string layer_name[] = {"SiliconSensor", "Metalconnection", "HDI", "Cooling", "Support", "Support_Gap", "Support2"};
    if (it.second == "SiliconSensor")
    {
      visatt->SetVisibility(true);
      visatt->SetColour(G4Colour(G4Colour::Green()));
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Metalconnection")
    {
      visatt->SetColour(105./255,105./255,94./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "HDI")
    {
      visatt->SetColour(85./255,85./255,74./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Support" || it.second == "Support2")
    {
      visatt->SetColour(160./255,160./255,200./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Support_Gap")
    {
      visatt->SetVisibility(false);
    }
    else if (it.second == "Cooling")
    {
      visatt->SetColour(G4Colour(0.823, 0.992, 0.980));
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
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
