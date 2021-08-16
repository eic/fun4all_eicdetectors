// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4HYBRIDHOMOGENEOUSCALORIMETERSUBSYSTEM_H
#define G4DETECTORS_PHG4HYBRIDHOMOGENEOUSCALORIMETERSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>  // for string

class PHCompositeNode;
class PHG4HybridHomogeneousCalorimeterDetector;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;

class PHG4HybridHomogeneousCalorimeterSubsystem : public PHG4DetectorSubsystem
{
 public:
  /** Constructor
   */
  PHG4HybridHomogeneousCalorimeterSubsystem(const std::string &name = "CRYSTAL_DEFAULT", const int layer = 0);

  /** Destructor
   */
  virtual ~PHG4HybridHomogeneousCalorimeterSubsystem();

  /**
     Creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
     Creates the stepping action and place it on the node tree, under "ACTIONS" node
     Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *) override;

  /** Event processing
   */
  int process_event(PHCompositeNode *) override;

  /** Accessors (reimplemented)
   */
  PHG4Detector *GetDetector(void) const override;
  PHG4SteppingAction *GetSteppingAction() const override { return m_SteppingAction; }
  PHG4DisplayAction *GetDisplayAction() const override { return m_DisplayAction; }

  /** Set mapping file for calorimeter towers
   */
  void SetTowerMappingFile(const std::string &filename);

  /** Select projective geometry for calorimeter
   */
  void SetProjectiveGeometry(const std::string &filename1, const std::string &filename2);

 private:
  //! set detector specific parameters and their defaults
  /*! called by PHG4DetectorSubsystem */
  void SetDefaultParameters() override;

  //  GeometryConfiguration current_geom_config_;

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4HybridHomogeneousCalorimeterDetector *m_Detector = nullptr;

  /** Stepping action
   */
  PHG4SteppingAction *m_SteppingAction = nullptr;
  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction = nullptr;
};

#endif
