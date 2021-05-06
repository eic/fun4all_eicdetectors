#include "G4LBLVtxDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>

#include <Geant4/G4Colour.hh>  // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4VisAttributes.hh>

using namespace std;

G4LBLVtxDisplayAction::G4LBLVtxDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

G4LBLVtxDisplayAction::~G4LBLVtxDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
  m_LogVolSet.clear();
}

void G4LBLVtxDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    string material_name(it->GetMaterial()->GetName());
    //    cout << "materialname: " << material_name << endl;
    if (material_name == "ITS_WATER")
    {
      VisAtt->SetColor(G4Colour::Blue());
    }
    else if (material_name == "ITS_PEEKCF30")
    {
      VisAtt->SetColor(G4Colour::White());
    }
    else if (material_name == "ITS_FGS003")
    {
      VisAtt->SetColor(G4Colour::Yellow());
    }
    else if (material_name == "ITS_CarbonFleece")
    {
      VisAtt->SetColor(G4Colour::Red());
    }
    else
    {
      VisAtt->SetColor(G4Colour::Cyan());
    }

    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(true);
    it->SetVisAttributes(VisAtt);
  }
  return;
}
