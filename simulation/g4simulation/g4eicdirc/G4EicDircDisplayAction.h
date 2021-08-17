// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4EICDIRCDISPLAYACTION_H
#define G4EICDIRCDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <string>  // for string
#include <map>

class G4Colour;
class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHParameters;

class G4EicDircDisplayAction : public PHG4DisplayAction
{
 public:
  G4EicDircDisplayAction(const std::string &name, PHParameters *parameters);

  virtual ~G4EicDircDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void SetMyVolume(G4LogicalVolume *vol) { m_MyVolume = vol; }
  void SetColor(const double red, const double green, const double blue, const double alpha = 1.);
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  PHParameters *m_Params;
  G4LogicalVolume *m_MyVolume;
  G4VisAttributes *m_VisAtt;
  G4Colour *m_Colour;
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
};

#endif  // G4DETECTORS_PHG4BLOCKDISPLAYACTION_H
