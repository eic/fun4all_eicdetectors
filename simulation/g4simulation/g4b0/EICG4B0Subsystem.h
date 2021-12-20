// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4B0SUBSYSTEM_H
#define EICG4B0SUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

class PHCompositeNode;
class PHG4Detector;
class EICG4B0Detector;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   *
   * The detector is constructed and registered via EICG4B0Detector
   *
   *
   * \see EICG4B0Detector
   * \see EICG4B0Subsystem
   *
   */
class EICG4B0Subsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  EICG4B0Subsystem(const std::string& name = "EICG4B0", const int layer = 0);

  //! destructor
  virtual ~EICG4B0Subsystem() {}

  /*!
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*) override;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector() const override;

  PHG4SteppingAction* GetSteppingAction() const override { return m_SteppingAction; }
  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

  bool CanBeMotherSubsystem() const override { return true; }

  void SaveAllHits(bool i = true) { m_SaveAllHitsFlag = i; }

  void SetTowerMappingFile(const std::string& filename);

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  EICG4B0Detector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

  bool m_SaveAllHitsFlag = false;

  std::string mappingfile_;
};

#endif  // EICG4B0SUBSYSTEM_H
