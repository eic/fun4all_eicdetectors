// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef G4DETECTORS_G4LBLVTXSUBSYSTEM_H
#define G4DETECTORS_G4LBLVTXSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>

class PHCompositeNode;
class G4LBLVtxDetector;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;

/*!
 * \brief G4LBLVtxSubsystem is a generic detector built from a GDML import
 */
class G4LBLVtxSubsystem : public PHG4DetectorSubsystem
{
 public:
  G4LBLVtxSubsystem(const std::string& name);
  virtual ~G4LBLVtxSubsystem();

  /*!
  creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*);

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const;
  PHG4SteppingAction* GetSteppingAction(void) const { return m_SteppingAction; }
  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  G4LBLVtxDetector* m_Detector;

  //! detector "stepping" action, executes after every G4 step
  /*! derives from PHG4SteppingAction */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;
};

#endif /* G4LBLVTXSUBSYSTEM_H_ */
