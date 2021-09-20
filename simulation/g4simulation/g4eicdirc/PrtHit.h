// -----------------------------------------
// PrtHit.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtHit_h
#define PrtHit_h 1

#include <Rtypes.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include "TVector3.h"

/*namespace bar_vectors
{
  std::vector<TVector3> vector_p_bar;
  std::vector<TVector3> vector_hit_pos_bar;

  std::vector<TVector3> get_p_bar() { return vector_p_bar; }
  std::vector<TVector3> get_pos_bar() { return vector_hit_pos_bar;}  

  }*/

class PrtHit : public PHG4Hitv1
{
 public:
  //Constructor
  PrtHit();

  ~PrtHit(){};

  // Accessors
  Int_t GetParticleId() { return fParticleId; }
  Int_t GetParentParticleId() { return fParentParticleId; }

  Int_t GetNreflectionsInPrizm() { return fNreflectionsInPrizm; }
  Long64_t GetPathInPrizm() { return fPathInPrizm; }

  TVector3 GetLocalPos() { return fLocalPos; }
  TVector3 GetGlobalPos() { return fGlobalPos; }
  TVector3 GetDigiPos() { return fDigiPos; }
  TVector3 GetMomentum() { return fMomentum; }
  TVector3 GetPosition() { return fPosition; }
  TVector3 GetMomentumAtBar() const { return fMomBar; };
  TVector3 GetPositionAtBar() const { return fPosBar; };

  Int_t GetMcpId() { return fMcpId; }
  Int_t GetPixelId() { return fPixelId; }
  Int_t GetChannel() { return fChannel; }
  Int_t GetMultiplicity() { return fMultiplicity; }
  Double_t GetLeadTime() { return fLeadTime; }
  Double_t GetTotTime() { return fTotTime; }

  // Mutators
  void SetParticleId(Int_t val) { fParticleId = val; }
  void SetParentParticleId(Int_t val) { fParentParticleId = val; }

  void SetNreflectionsInPrizm(Int_t val) { fNreflectionsInPrizm = val; }
  void SetPathInPrizm(Long64_t val) { fPathInPrizm = val; }

  void SetLocalPos(TVector3 val) { fLocalPos = val; }
  void SetGlobalPos(TVector3 val) { fGlobalPos = val; }
  void SetDigiPos(TVector3 val) { fDigiPos = val; }
  void SetMomentum(TVector3 val) { fMomentum = val; }
  void SetPosition(TVector3 val) { fPosition = val; }
  void SetMomentumAtBar(TVector3 val) { fMomBar = val; }
  void SetPositionAtBar(TVector3 val) { fPosBar = val; }

  void SetMcpId(Int_t val) { fMcpId = val; }
  void SetPixelId(Int_t val) { fPixelId = val; }
  void SetChannel(Int_t val) { fChannel = val; }
  void SetMultiplicity(Int_t val) { fMultiplicity = val; }
  void SetLeadTime(Double_t val) { fLeadTime = val; }
  void SetTotTime(Double_t val) { fTotTime = val; }

 protected:
  Int_t fParticleId;
  Int_t fParentParticleId;
  Int_t fNreflectionsInPrizm;
  Long64_t fPathInPrizm;
  TVector3 fLocalPos;
  TVector3 fGlobalPos;
  TVector3 fDigiPos;
  TVector3 fMomentum;
  TVector3 fPosition;
  TVector3 fMomBar;
  TVector3 fPosBar;

  Int_t fMcpId;
  Int_t fPixelId;
  Int_t fChannel;
  Int_t fMultiplicity;
  Double_t fLeadTime;
  Double_t fTotTime;

  ClassDef(PrtHit, 2)
};

#endif
