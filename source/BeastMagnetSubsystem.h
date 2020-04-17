// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BEASTMAGNETSUBSYSTEM_H
#define BEASTMAGNETSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

class BeastMagnetDetector;
class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   *
   * The detector is constructed and registered via BeastMagnetDetector
   *
   *
   * \see BeastMagnetDetector
   * \see BeastMagnetSubsystem
   *
   */
class BeastMagnetSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  BeastMagnetSubsystem(const std::string& name = "BeastMagnet");

  //! destructor
  virtual ~BeastMagnetSubsystem();

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
  PHG4DisplayAction* GetDisplayAction() const  override { return m_DisplayAction; }

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  BeastMagnetDetector  *m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;
};

#endif // BEASTMAGNETSUBSYSTEM_H
