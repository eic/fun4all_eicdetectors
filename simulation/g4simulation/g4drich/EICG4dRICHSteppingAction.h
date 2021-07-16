#ifndef DRICHSTEPPINGACTION_H
#define DRICHSTEPPINGACTION_H

#include <G4StepPoint.hh>
#include <G4String.hh>
#include <G4Track.hh>
#include <g4main/PHG4SteppingAction.h>

class EICG4dRICHDetector;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class EICG4dRICHHit;
class PHG4HitContainer;
class PHParameters;

class EICG4dRICHSteppingAction : public PHG4SteppingAction 
{
  public:
    //! constructor
    EICG4dRICHSteppingAction(EICG4dRICHDetector *, const PHParameters *parameters);

    //! destructor
    virtual ~EICG4dRICHSteppingAction();

    //! stepping action
    virtual bool UserSteppingAction(const G4Step *, bool);

    //! reimplemented from base class
    virtual void SetInterfacePointers(PHCompositeNode *);

  private:
    //! method to initialize a new hit, resetting some things, such
    //  as energy deposition accumulators
    void InitHit(const G4StepPoint *prePoint_, const G4Track *aTrack_,
                 bool resetAccumulators);

    //! pointer to the detector
    EICG4dRICHDetector *m_Detector;
    const PHParameters *m_Params;
    //! pointer to hit container
    PHG4HitContainer *m_HitContainer;
    EICG4dRICHHit *m_Hit;
    PHG4HitContainer *m_SaveHitContainer;
    G4VPhysicalVolume *m_SaveVolPre;
    G4VPhysicalVolume *m_SaveVolPost;

    int m_SaveTrackId;
    int m_SavePreStepStatus;
    int m_SavePostStepStatus;
    int m_ActiveFlag;
    double m_EdepSum;
    double m_EionSum;

    // hit type classifiers
    int hitType;
    enum hitTypes 
    {
      hEntrance, /* vessel entrance */
      hExit,     /* vessel exit */
      hPSST,     /* photosensor hit */
      hIgnore,   /* none of the above */
      nHitTypes
    };
    G4String hitTypeStr[nHitTypes];
    // hit subtype classifiers
    int hitSubtype;
    enum hitSubtypes 
    {
      /* entrance hits                    */
      entPrimary,   /* primary, thrown from generator */
      entSecondary, /* secondary, byproduct of thrown particle */
      entPostStep,  /* incident particle from PostStepDoItVector */
      /* exit hits                        */
      exPrimary,   /* primary track exit */
      exSecondary, /* secondary track exit (not primary) */
      /* photosensor hits                 */
      psOptical, /* opticalphoton hit */
      psGamma,   /* non-optical photon hit */
      psOther,   /* non-photon hit */
      /* unknown hit                      */
      subtypeUnknown,
      nHitSubtypes
    };
    G4String hitSubtypeStr[nHitSubtypes];
};

#endif // DRICHSTEPPINGACTION_H
