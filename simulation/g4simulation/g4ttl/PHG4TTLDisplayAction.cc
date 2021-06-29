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

PHG4TTLDisplayAction::PHG4TTLDisplayAction(const std::string &name, bool detailed)
  : PHG4DisplayAction(name)
{
  showdetails = detailed;
  if (!detailed) std::cout << "PHG4TTLDisplayAction::disabled detailed view of towers" << std::endl;
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
      if (!showdetails) visatt->SetVisibility(false);
    }
    else if (it.second == "DetectorBox")
    {
      // visatt->SetVisibility(false);
      
      if (!showdetails){
        visatt->SetColour(G4Colour::Black());
        visatt->SetForceWireframe(true);
      } else {
        visatt->SetColour(34./255, 139./255, 34./255, 1. );
      }
      visatt->SetForceLineSegmentsPerCircle(50);
      
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
