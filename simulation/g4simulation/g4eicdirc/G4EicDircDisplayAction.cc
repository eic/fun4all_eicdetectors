#include "G4EicDircDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Utils.h>

#include <phparameter/PHParameters.h>

#include <Geant4/G4Colour.hh>  // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>  // for isfinite

G4EicDircDisplayAction::G4EicDircDisplayAction(const std::string &name, PHParameters *pars)
  : PHG4DisplayAction(name)
  , m_Params(pars)
  , m_MyVolume(nullptr)
  , m_VisAtt(nullptr)
  , m_Colour(nullptr)
{
}

G4EicDircDisplayAction::~G4EicDircDisplayAction()
{
  delete m_VisAtt;
  delete m_Colour;
}

void G4EicDircDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing
  for (auto it : m_LogicalVolumeMap)
  {
    m_MyVolume = it.first;
    // if (m_MyVolume->GetVisAttributes())
    // {
    //   m_VisAtt = new G4VisAttributes();
    //   m_VisAtt->SetForceWireframe(false);
    //   m_VisAtt->SetForceSolid(true);
    //   continue;
    // }
    m_VisAtt = new G4VisAttributes();
    if (m_Params->get_int_param("blackhole"))
    {
      PHG4Utils::SetColour(m_VisAtt, "BlackHole");
      m_VisAtt->SetVisibility(false);
      m_VisAtt->SetForceSolid(false);
    }
    else
    {
      PHG4Utils::SetColour(m_VisAtt, m_Params->get_string_param("material"));
      m_VisAtt->SetVisibility(true);
      m_VisAtt->SetForceWireframe(false);
      m_VisAtt->SetForceSolid(true);
    }
    if (m_Colour)
    {
      m_VisAtt->SetColour(m_Colour->GetRed(),
                          m_Colour->GetGreen(),
                          m_Colour->GetBlue(),
                          m_Colour->GetAlpha());
      m_VisAtt->SetForceWireframe(false);
      m_VisAtt->SetVisibility(true);
      m_VisAtt->SetForceSolid(true);
    }
    m_MyVolume->SetVisAttributes(m_VisAtt);
  }
  return;
}

void G4EicDircDisplayAction::SetColor(const double red, const double green, const double blue, const double alpha)
{
  if (std::isfinite(red) && std::isfinite(green) && std::isfinite(blue) && std::isfinite(alpha))
  {
    m_Colour = new G4Colour(red, green, blue, alpha);
  }
  return;
}
