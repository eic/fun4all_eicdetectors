#ifndef G4EICDIRCSTACKINGACTION_H
#define G4EICDIRCSTACKINGACTION_H

#include <g4main/PHG4StackingAction.h>

#include <Geant4/G4UserStackingAction.hh>
#include <Geant4/globals.hh>

#include <gsl/gsl_rng.h>

class G4EicDircDetector;
class TGraph;

class G4EicDircStackingAction : public PHG4StackingAction
{
 public:
  G4EicDircStackingAction(G4EicDircDetector* det);
  ~G4EicDircStackingAction() override;

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack) override;
  void PrepareNewEvent() override;


 private:
  gsl_rng* m_RandomGenerator = nullptr;
  G4EicDircDetector* m_Detector = nullptr;
  TGraph* fDetEff[2] = {nullptr};
  int fCerenkovCounter = 0;
  int fScintillationCounter = 0;
};

#endif
