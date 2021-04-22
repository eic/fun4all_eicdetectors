// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ALLSILICONTRACKERSUBSYSTEM_H
#define ALLSILICONTRACKERSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <set>
#include <string>

class AllSiliconTrackerDetector;
class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   *
   * The detector is constructed and registered via AllSiliconTrackerDetector
   *
   *
   * \see AllSiliconTrackerDetector
   * \see AllSiliconTrackerSubsystem
   *
   */
class AllSiliconTrackerSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  AllSiliconTrackerSubsystem(const std::string& name = "AllSiliconTracker");

  //! destructor
  virtual ~AllSiliconTrackerSubsystem();

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
  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;
  void AddAssemblyVolume(const std::string& avol);
  std::pair<std::set<std::string>::const_iterator, std::set<std::string>::const_iterator> assembly_iters() { return std::make_pair(m_AssemblyVolumeSet.begin(), m_AssemblyVolumeSet.end()); }

  void AddLogicalVolume(const std::string& name);
  std::pair<std::set<std::string>::const_iterator, std::set<std::string>::const_iterator> logvol_iters() { return std::make_pair(m_LogVolumeSet.begin(), m_LogVolumeSet.end()); }

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  AllSiliconTrackerDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  std::set<std::string> m_AssemblyVolumeSet;
  std::set<std::string> m_LogVolumeSet;
};

#endif  // ALLSILICONTRACKERSUBSYSTEM_H
