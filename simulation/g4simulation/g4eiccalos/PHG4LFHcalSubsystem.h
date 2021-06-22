// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4LFHCALSUBSYSTEM_H
#define G4DETECTORS_PHG4LFHCALSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>  // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4LFHcalDetector;
class PHG4SteppingAction;

class PHG4LFHcalSubsystem : public PHG4DetectorSubsystem
{
 public:
  /** Constructor
   */
  PHG4LFHcalSubsystem(const std::string &name = "LF_HCAL", const int layer = 0);

  /** Destructor
   */
  virtual ~PHG4LFHcalSubsystem();

  /**
     Creates the m_Detector object
     Creates the stepping action
     Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *);

  /** Event processing
   */
  int process_event(PHCompositeNode *);

  /** Accessors (reimplemented)
   */
  PHG4Detector *GetDetector() const;
  PHG4SteppingAction *GetSteppingAction() const { return m_SteppingAction; }
  PHG4DisplayAction *GetDisplayAction() const { return m_DisplayAction; }

  /** Set mapping file for calorimeter towers
   */
  void SetTowerMappingFile(const std::string &filename);
  /** Set layers per tower segment by hand, has to be set in addition to mapping file
   */
  void SetLayerPerTowerSegment(int layerPerTowerSeg) {set_int_param("nlayerspertowerseg", layerPerTowerSeg); };
  
 private:
  void SetDefaultParameters();

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4LFHcalDetector *m_Detector = nullptr;

  /** Stepping action
   */
  PHG4SteppingAction *m_SteppingAction = nullptr;
  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction = nullptr;
};

#endif
