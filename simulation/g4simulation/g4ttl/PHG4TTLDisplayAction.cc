#include "PHG4TTLDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Utils.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

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
    else if (it.second == "TTLLayers")
    {
      // visatt->SetForceWireframe(true);
      PHG4Utils::SetColour(visatt, it.first->GetMaterial()->GetName());
      // if(it.first->GetMaterial()->GetName()=="CarbonFiberSupport") visatt->SetColour(G4Colour(21./255, 27./255, 31./255));
      if (it.first->GetMaterial()->GetName() == "CarbonFiberSupport"){
        visatt->SetVisibility(false);
        visatt->SetColour(G4Colour(4 * 21. / 255, 4 * 27. / 255, 4 * 31. / 255));
      }
      if (it.first->GetMaterial()->GetName() == "G4_Al"){
        visatt->SetVisibility(false);
        visatt->SetColour(G4Colour(132. / 255, 135. / 255, 137. / 255));
      }
      if (it.first->GetMaterial()->GetName() == "G4_GRAPHITE"){
        visatt->SetVisibility(false);
        visatt->SetColour(G4Colour(2 * 21. / 255, 2 * 27. / 255, 2 * 31. / 255));
      }
      if (it.first->GetMaterial()->GetName() == "AluminiumNitrate"){
        visatt->SetVisibility(false);
        visatt->SetColour(G4Colour(0.8 * 138, 0.8 * 115, 0.8 * 115));
      }
      if (it.first->GetMaterial()->GetName() == "G4_PLEXIGLASS"){
        visatt->SetVisibility(false);
        visatt->SetColour(G4Colour(68, 131, 157));
      }
    }
    else if (it.second == "SHLayers")
    {
      // visatt->SetForceWireframe(true);
      PHG4Utils::SetColour(visatt, it.first->GetMaterial()->GetName());
      if (it.first->GetMaterial()->GetName() == "G4_GRAPHITE") visatt->SetColour(G4Colour(2 * 21. / 255, 2 * 27. / 255, 2 * 31. / 255));
      if (it.first->GetMaterial()->GetName() == "G4_POLYSTYRENE") visatt->SetColour(G4Colour(193, 89, 0));
    }
    else if (it.second == "Support")
    {
      // visatt->SetColour(G4Colour::White());
      visatt->SetColour(G4Colour::Gray());
      // visatt->SetColour(G4Colour(4*21./255, 4*27./255, 4*31./255,0.8));
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "StripBox")
    {
      visatt->SetColour(G4Colour::Blue());
      // visatt->SetForceSolid(true);
      visatt->SetForceWireframe(true);
      visatt->SetVisibility(false);
    }
    else if (it.second == "Cooling_tube")
    {
      visatt->SetColour(G4Colour(0.7, 0.7, 0.7));
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
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
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "SensorStack" || it.second == "SensorLadder" || it.second == "SensorAndReadoutLadder")
    {
      visatt->SetVisibility(false);
      visatt->SetColour(G4Colour(G4Colour::Green()));
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "SHLadder" || it.second == "SHStack")
    {
      visatt->SetVisibility(false);
      visatt->SetColour(G4Colour(G4Colour::Green()));
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "ModuleEnvelope")
    {
      // visatt->SetVisibility(false);
      visatt->SetColour(G4Colour::Red());
      // visatt->SetForceSolid(true);
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "FullEnvelope")
    {
      visatt->SetVisibility(false);
      visatt->SetColour(G4Colour::Red());
      // visatt->SetForceSolid(true);
      visatt->SetForceWireframe(true);
    }

    else if (it.second == "Water_cooling")
    {
      visatt->SetColour(G4Colour(0.823, 0.992, 0.980));
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "Cooling_Support")
    {
      // visatt->SetColour(G4Colour(21./255, 27./255, 31./255,0.5));
      visatt->SetColour(G4Colour(4 * 21. / 255, 4 * 27. / 255, 4 * 31. / 255));
      visatt->SetForceSolid(true);
      // visatt->SetForceWireframe(true);
    }
    else if (it.second == "StripCoolingSupportBox")
    {
      visatt->SetColour(G4Colour::Green());
      // visatt->SetForceSolid(true);
      visatt->SetForceWireframe(true);
      visatt->SetVisibility(false);
    }
    else if (it.second == "Module_Mother")
    {
      visatt->SetVisibility(false);
    }
    else if (it.second == "DetectorBox")
    {
      // visatt->SetVisibility(false);

      // if (showdetails){
      // visatt->SetColour(G4Colour::Black());
      visatt->SetColour(G4Colour::Red());
      visatt->SetForceWireframe(true);
      visatt->SetVisibility(false);
      visatt->SetForceLineSegmentsPerCircle(50);
      // } else {
      //   visatt->SetColour(34./255, 139./255, 34./255, 1. );
      //   visatt->SetForceSolid(true);
      //   visatt->SetVisibility(true);
      // }
    }
    else if (it.second == "DetectorBoxFwd")
    {
      // if (showdetails){
      visatt->SetColour(G4Colour::Red());
      visatt->SetForceWireframe(true);
      visatt->SetVisibility(true);
      visatt->SetForceLineSegmentsPerCircle(50);
      // } else {
      //   visatt->SetColour(34./255, 139./255, 34./255, 1. );
      //   visatt->SetForceSolid(true);
      //   visatt->SetVisibility(true);
      // }
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
