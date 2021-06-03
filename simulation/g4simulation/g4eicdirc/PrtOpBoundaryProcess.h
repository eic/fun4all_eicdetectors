#ifndef PrtOpBoundaryProcess_h
#define PrtOpBoundaryProcess_h

#include <Geant4/globals.hh>
#include <Geant4/G4OpBoundaryProcess.hh>

class PrtOpBoundaryProcess : public G4OpBoundaryProcess
{
public:
  PrtOpBoundaryProcess();
  ~PrtOpBoundaryProcess() override {};

public:
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) override;

private:
  int fLensId;
};


#endif /*PrtOpBoundaryProcess_h*/
