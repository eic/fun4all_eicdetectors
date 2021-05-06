// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4LBLVTXDISPLAYACTION_H
#define G4LBLVTXDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <set>
#include <string>  // for string
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class G4LBLVtxDisplayAction : public PHG4DisplayAction
{
 public:
  G4LBLVtxDisplayAction(const std::string &name);

  virtual ~G4LBLVtxDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void AddLogVolume(G4LogicalVolume *vol) { m_LogVolSet.insert(vol); }

 private:
  std::vector<G4VisAttributes *> m_VisAttVec;
  std::set<G4LogicalVolume *> m_LogVolSet;
};

#endif  // G4LBLVTXDISPLAYACTION_H
