// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDHCALDISPLAYACTION_H
#define G4DETECTORS_PHG4FORWARDHCALDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class PHG4ForwardHcalDisplayAction : public PHG4DisplayAction
{
 public:
  explicit PHG4ForwardHcalDisplayAction(const std::string &name);
  explicit PHG4ForwardHcalDisplayAction(const std::string &name, bool detailed);

  virtual ~PHG4ForwardHcalDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
  bool showdetails = false;
};

#endif  // G4DETECTORS_PHG4FORWARDHCALDISPLAYACTION_H
