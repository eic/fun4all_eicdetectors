// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4EICDIRCSUBSYSTEM_H
#define G4EICDIRCSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>

class G4EicDircDetector;
class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;

/**
   * \brief example Fun4All module
   *
   * The detector is constructed and registered via G4EicDircDetector
   *
   *
   * \see G4EicDircDetector
   * \see G4EicDircSubsystem
   *
   */
class G4EicDircSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  G4EicDircSubsystem(const std::string& name = "G4EicDirc");

  //! destructor
  virtual ~G4EicDircSubsystem();

  /*!
  creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  virtual int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  virtual int process_event(PHCompositeNode*);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector() const;

  virtual PHG4SteppingAction* GetSteppingAction() const { return m_SteppingAction; }
  //! Print info (from SubsysReco)
  virtual void Print(const std::string& what = "ALL") const;

  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

 private:
  // \brief Set default parameter values
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  G4EicDircDetector* m_Detector;

  PHG4SteppingAction* m_SteppingAction;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;
  //! Color setting if we want to override the default

  std::string m_HitNodeName;
  std::string m_AbsorberNodeName;
  std::string m_SupportNodeName;

};

#endif  // G4EICDIRCSUBSYSTEM_H
