#include "AllSiliconTrackerDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Utils.h>

#include <Geant4/G4Colour.hh>  // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>  // for isfinite

using namespace std;

AllSiliconTrackerDisplayAction::AllSiliconTrackerDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
  , m_MyVolume(nullptr)
{
}

AllSiliconTrackerDisplayAction::~AllSiliconTrackerDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
  m_LogVolSet.clear();
}

void AllSiliconTrackerDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing
  for (auto &it : m_LogVolSet)
  {
    if (it->GetVisAttributes())
    {
      continue;
    }

    //    cout << "mat: " << it->GetMaterial()->GetName() << endl;
    G4VisAttributes *VisAtt = new G4VisAttributes();
    m_VisAttVec.push_back(VisAtt);
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(true);
    if (it->GetMaterial()->GetName().find("air") != string::npos)
    {
      VisAtt->SetVisibility(false);
      VisAtt->SetForceSolid(false);
    }
    else if (it->GetMaterial()->GetName().find("CarbonFiber") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Grey());
    }
    else if (it->GetMaterial()->GetName().find("Kapton") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Green());
    }
    else if (it->GetMaterial()->GetName().find("water") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Blue());
    }
    else if (it->GetMaterial()->GetName().find("silicon") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Cyan());
    }
    else if (it->GetMaterial()->GetName().find("aluminum") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Yellow());
    }
    else if (it->GetMaterial()->GetName().find("beryllium") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Red());
    }
    else
    {
      VisAtt->SetColour(G4Colour::Red());
    }

    it->SetVisAttributes(VisAtt);
  }
  return;
}
