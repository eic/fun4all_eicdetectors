// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4ZDCSUBSYSTEM_H
#define EICG4ZDCSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

class PHCompositeNode;
class PHG4Detector;
class EICG4ZDCDetector;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   *
   * The detector is constructed and registered via EICG4ZDCDetector
   *
   *
   * \see EICG4ZDCDetector
   * \see EICG4ZDCSubsystem
   *
   */
class EICG4ZDCSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  EICG4ZDCSubsystem(const std::string& name = "EICG4ZDC");

  //! destructor
  virtual ~EICG4ZDCSubsystem() {}

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

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  EICG4ZDCDetector  *m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *m_SteppingAction;
};

#endif // EICG4ZDCSUBSYSTEM_H
