// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ALLSI_AL_SUPPORT_DISPLAYACTION_H
#define ALLSI_AL_SUPPORT_DISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <set>
#include <string>  // for string
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class AllSi_Al_support_DisplayAction : public PHG4DisplayAction
{
 public:
  AllSi_Al_support_DisplayAction(const std::string &name);

  virtual ~AllSi_Al_support_DisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void SetMyVolume(G4LogicalVolume *vol) { m_MyVolume = vol; }
  void AddLogicalVolume(G4LogicalVolume *vol) { m_LogVolSet.insert(vol); }

 private:
  G4LogicalVolume *m_MyVolume;
  std::vector<G4VisAttributes *> m_VisAttVec;
  std::set<G4LogicalVolume *> m_LogVolSet;
};

#endif  // ALLSI_AL_SUPPORT_DISPLAYACTION_H
