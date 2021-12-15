#include "PHG4LFHcalDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>                     // for pair

using namespace std;

PHG4LFHcalDisplayAction::PHG4LFHcalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4LFHcalDisplayAction::PHG4LFHcalDisplayAction(const std::string &name, bool detailed)
  : PHG4DisplayAction(name)
{
  showdetails = detailed;
  if (!detailed) std::cout << "PHG4LFHcalDisplayAction::disabled detailed view of towers" << std::endl;
}

PHG4LFHcalDisplayAction::~PHG4LFHcalDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4LFHcalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing
  for (auto it : m_LogicalVolumeMap)
  {
    G4LogicalVolume *logvol = it.first;
    if (logvol->GetVisAttributes())
    {
      continue;
    }
    // std::cout << showdetails << "\t" << it.second ;
    G4VisAttributes *visatt = new G4VisAttributes();
    visatt->SetForceSolid(true);
    m_VisAttVec.push_back(visatt);  // for later deletion
    if (it.second == "Invisible")
    {
        visatt->SetVisibility(false);
    }
    else if (it.second == "Wireframe")
    {
        visatt->SetColour(G4Colour::Red());
        visatt->SetForceWireframe(true);
        visatt->SetVisibility(true);
    }
    else if (it.second == "Absorber_W")
    {
      if (showdetails){
        visatt->SetColour(G4Colour::Black());
        visatt->SetVisibility(true);
        // visatt->SetForceWireframe(true);
        // std::cout << " is visible" ;
      } else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "Absorber")
    {
      if (showdetails){
        visatt->SetColour(4 * 21. / 255, 4 * 27. / 255, 4 * 31. / 255);
        visatt->SetVisibility(true);
        // visatt->SetForceWireframe(true);
        // std::cout << " is visible" ;
      } else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "Spacer")
    {
      if (showdetails){
        visatt->SetColour(220. / 255, 220. / 255, 220. / 255);
        visatt->SetVisibility(true);
        // visatt->SetForceWireframe(true);
        // std::cout << " is visible" ;
      } else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "LFHcalEnvelope")
    {
      visatt->SetColour(G4Colour::Magenta());
      visatt->SetForceSolid(false);
      visatt->SetVisibility(false);
    }
    else if (it.second == "Scintillator")
    {
      if (showdetails){
        // visatt->SetColour(G4Colour::White());
        visatt->SetColour(127./ 255,255./ 255,212./ 255, 0.5);
        visatt->SetVisibility(true);
        // visatt->SetForceWireframe(true);
        // std::cout << " is visible" ;
      } else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "WLSfiber")
    {
      if (showdetails){
        // visatt->SetColour(G4Colour::Blue());
        visatt->SetColour(20./ 255,10./ 255,200./ 255, 1.0);
        // visatt->SetColour(152./ 255,10./ 255,10./ 255, 1.0);
        // visatt->SetColour(152./ 255,251./ 255,152./ 255, 0.9);
        visatt->SetVisibility(true);
        // std::cout << " is visible" ;
      } else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "Frame")
    {
      if (showdetails){
        visatt->SetColour(G4Colour::Yellow());
        // visatt->SetColour(152./ 255,10./ 255,10./ 255, 0.9);
        // visatt->SetColour(152./ 255,251./ 255,152./ 255, 0.9);
        visatt->SetVisibility(true);
        visatt->SetForceWireframe(true);
        // std::cout << " is visible" ;
      } else 
        visatt->SetVisibility(false);
    }
    else if (it.second == "SingleTower_W")
    {
      if (showdetails)
        visatt->SetVisibility(false);
      else {
        visatt->SetColour(G4Colour::Black());
        visatt->SetVisibility(true);
        // std::cout << " is visible" ;
      }
    }
    else if (it.second == "SingleTower")
    {
      if (showdetails)
        visatt->SetVisibility(false);
      else {
        visatt->SetColour(G4Colour::Blue());
        visatt->SetVisibility(true);
        // std::cout << " is visible" ;
      }
    }
    else
    {
      cout << "unknown logical volume " << it.second << endl;
      gSystem->Exit(1);
    }
    // std::cout << std::endl ;
    logvol->SetVisAttributes(visatt);
  }
  return;
}
