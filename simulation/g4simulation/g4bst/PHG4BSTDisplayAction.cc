#include "PHG4BSTDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>                     // for pair

using namespace std;

PHG4BSTDisplayAction::PHG4BSTDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4BSTDisplayAction::~PHG4BSTDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4BSTDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    if (it.second == "Absorber")
    {
      visatt->SetColour(G4Colour::Gray());
      visatt->SetVisibility(true);
    }
    else if (it.second == "InnerBarrel")
    {
      // visatt->SetColour(212./255,175./255,55./255);
      visatt->SetColour(212./255,175./255,55./255,0.6);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "BackingMaterial")
    {
      visatt->SetColour(212./255,50./255,55./255);
      visatt->SetVisibility(true);
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "CopperWire")
    {
      visatt->SetColour(G4Color::Red());
      // visatt->SetColour(250./255,10./255,10./255, 0.7);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "OuterBarrel")
    {
      visatt->SetColour(207./255,181./255,59./255,0.6);
      // visatt->SetColour(207./255,181./255,59./255);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Support")
    {
      visatt->SetColour(62./255,62./255,64./255);
      visatt->SetVisibility(true);
    }
    else if (it.second == "FoamEndWheel")
    {
      visatt->SetColour(0.6*85./255,0.6*85./255,0.6*74./255, 1.0);
      visatt->SetVisibility(true);
    }
    else if (it.second == "Foam")
    {
      visatt->SetColour(85./255,85./255,74./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "CarbonPlate")
    {
      visatt->SetColour(105./255,105./255,94./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "CarbonFleece")
    {
      visatt->SetColour(85./255,85./255,74./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "GraphiteFoil")
    {
      visatt->SetColour(160./255,160./255,200./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Glue")
    {
      visatt->SetColour(G4Colour(G4Colour::White()));
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Kapton")
    {
      visatt->SetColour(G4Colour(255. / 255, 165. / 255, 0. / 255));
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "CircuitBoard")
    {
      visatt->SetColour(G4Colour(132. / 255, 135. / 255, 137. / 255));  // Al color
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "CarbonFoam")
    {
      visatt->SetColour(105./255,105./255,94./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Cooling_tube")
    {
      visatt->SetColour(G4Colour(0.7, 0.7, 0.7));
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
    }

    else if (it.second == "Water_cooling")
    {
      visatt->SetColour(G4Colour(0.823, 0.992, 0.980));
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Cooling_Support" || it.second == "Carbon_Support")
    {
      // visatt->SetColour(G4Colour(21./255, 27./255, 31./255,0.5));
      visatt->SetColour(G4Colour(4 * 21. / 255, 4 * 27. / 255, 4 * 31. / 255));
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
      visatt->SetVisibility(false);
    }
    else if (it.second == "StripCoolingSupportBox")
    {
      visatt->SetColour(G4Colour::Green());
      // visatt->SetForceSolid(true);
      visatt->SetForceWireframe(true);
      visatt->SetVisibility(false);
    }
    // NEW TTL:
    else if (it.second == "CoolingPlate")
    {
      visatt->SetColour(G4Colour(132. / 255, 135. / 255, 137. / 255));  // Al color
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Epoxy")
    {
      visatt->SetColour(G4Colour(G4Colour::White()));
      // visatt->SetForceSolid(true);
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "Sensor")
    {
      visatt->SetVisibility(true);
      visatt->SetColour(G4Colour(G4Colour::Green()));
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "CarbonHoneycomb")
    {
      visatt->SetColour(75./255,75./255,64./255, 1.0);
      visatt->SetVisibility(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "CShell")
    {
      visatt->SetColour(G4Colour::Blue());
      // visatt->SetColour(G4Colour::Black());
      visatt->SetVisibility(true);
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "Invisible")
    {
      visatt->SetVisibility(false);
    }
    else if (it.second == "FbstEnvelope")
    {
      visatt->SetVisibility(true);
      visatt->SetForceWireframe(true);
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "TestScint")
    {
      visatt->SetColour(G4Colour::White());
    }
    else if (it.second == "TestAbsorb")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "Scintillator")
    {
      visatt->SetColour(G4Colour::White());
      visatt->SetVisibility(false);
    }
    else if (it.second == "Cherenkov")
    {
      visatt->SetColour(G4Colour::Yellow());
      visatt->SetVisibility(false);
    }
    else if (it.second == "SingleTowerAbsorber")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "StripBox")
    {
      visatt->SetColour(G4Colour::Blue());
      // visatt->SetForceSolid(true);
      visatt->SetForceWireframe(true);
      visatt->SetVisibility(true);
    }
    else if (it.second == "ASIC")
    {
      visatt->SetVisibility(true);
      visatt->SetColour(G4Colour(G4Colour::Red()));
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "SingleTowerCherenkovAbsorber")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "SingleTowerCherenkov")
    {
      visatt->SetColour(G4Colour::Cyan());
      visatt->SetVisibility(false);
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
