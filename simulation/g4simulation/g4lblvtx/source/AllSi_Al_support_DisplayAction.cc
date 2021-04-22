#include "AllSi_Al_support_DisplayAction.h"

#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Utils.h>

#include <Geant4/G4Colour.hh>  // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>  // for isfinite

using namespace std;

AllSi_Al_support_DisplayAction::AllSi_Al_support_DisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
  , m_MyVolume(nullptr)
{
}

AllSi_Al_support_DisplayAction::~AllSi_Al_support_DisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
  m_LogVolSet.clear();
}

void AllSi_Al_support_DisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    else if (it->GetMaterial()->GetName().find("aluminum") != string::npos)
    {
      VisAtt->SetColour(G4Colour::Yellow());
    }
    else
    {
      VisAtt->SetColour(G4Colour::Red());
    }

    it->SetVisAttributes(VisAtt);
  }
  return;
}
