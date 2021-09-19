// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4EICDIRCSTEPPINGACTION_H
#define G4EICDIRCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>
#include <vector>
#include <Rtypes.h>
#include <TVector3.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class G4EicDircDetector;
class PHG4Hit;
class PHG4Hitv1;
class PHG4HitContainer;
class PHParameters;
class PrtHit;
class G4Track;

class G4EicDircSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  G4EicDircSteppingAction(G4EicDircDetector*, const PHParameters* parameters);

  //! destructor
  ~G4EicDircSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

  void SetHitNodeName(const std::string& nam) { m_HitNodeName = nam; }
  void SetAbsorberNodeName(const std::string& nam) { m_AbsorberNodeName = nam; }
  void SetSupportNodeName(const std::string& nam) { m_SupportNodeName = nam; }

  
  std::vector<Int_t> vector_nid;
  std::vector<Int_t> vector_trackid;

  //std::vector<Int_t> vector_bar_hit_trackid;
  //std::vector<TVector3> vector_p_bar;
  //std::vector<TVector3> vector_hit_pos_bar;

 private:
  //! pointer to the detector
  G4EicDircDetector* m_Detector = nullptr;
  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer = nullptr;
  PHG4HitContainer* m_AbsorberHitContainer = nullptr;
  PHG4HitContainer* m_SupportHitContainer = nullptr;
  PHG4Hit* m_Hit = nullptr;
//  PrtHit* m_Hit = nullptr;
  PHG4HitContainer* m_SaveHitContainer = nullptr;

  G4VPhysicalVolume* m_SaveVolPre = nullptr;
  G4VPhysicalVolume* m_SaveVolPost = nullptr;
  int m_SaveTrackId = -1;
  int m_SavePreStepStatus = -1;
  int m_SavePostStepStatus = -1;
  int m_ActiveFlag = 0;
  int m_BlackHoleFlag = 0;
  double m_EdepSum = 0.;
  double m_EionSum = 0.;

  std::string m_HitNodeName;
  std::string m_AbsorberNodeName;
  std::string m_SupportNodeName;
};

#endif  // G4EICDIRCSTEPPINGACTION_H
