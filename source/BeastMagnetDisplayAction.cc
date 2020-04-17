#include "BeastMagnetDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Utils.h>

#include <Geant4/G4Colour.hh>  // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>  // for isfinite

using namespace std;

BeastMagnetDisplayAction::BeastMagnetDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
  , m_MyVolume(nullptr)
{
}

BeastMagnetDisplayAction::~BeastMagnetDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
  m_LogVolSet.clear();
}

void BeastMagnetDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing
  for (auto &it : m_LogVolSet)
  {
    if (it->GetVisAttributes())
    {
      continue;
    }
    G4VisAttributes *VisAtt = new G4VisAttributes();
    m_VisAttVec.push_back(VisAtt);
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(true);
    if (it->GetName().find("Coil") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Cyan());
    }
    else if (it->GetName().find("CryostatHe") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Blue());
    }
    else if (it->GetName().find("CryostatAl") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Green());
    }
    else
    {
      VisAtt->SetColour(G4Colour::Red());
    }

    it->SetVisAttributes(VisAtt);
  }
  return;
}
