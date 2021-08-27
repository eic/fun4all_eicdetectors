#ifndef G4EICDIRCOPBOUNDARYPROCESS_H
#define G4EICDIRCOPBOUNDARYPROCESS_H

#include <Geant4/G4OpBoundaryProcess.hh>
#include <Geant4/G4OpticalPhoton.hh>
#include <Geant4/G4VParticleChange.hh>

class G4EicDircOpBoundaryProcess : public G4OpBoundaryProcess
{
public:
  explicit G4EicDircOpBoundaryProcess(const G4String& processName = "G4EicDircOpBoundary",
				    G4ProcessType type = fOptical) : G4OpBoundaryProcess()
{}

  virtual ~G4EicDircOpBoundaryProcess() override {}

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) override;

private:

};

inline
G4bool G4EicDircOpBoundaryProcess::IsApplicable(const G4ParticleDefinition&
                                                       aParticleType)
{
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

#endif
