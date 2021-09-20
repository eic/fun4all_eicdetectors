// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERSTRIPSUBSYSTEM_H
#define G4DETECTORS_PHG4CYLINDERSTRIPSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#if !defined(__CINT__) || defined(__CLING__)
#include <array>   // for array
#endif

#include <string>  // for string

class PHCompositeNode;
class PHG4CylinderStripDetector;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;

class PHG4CylinderStripSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4CylinderStripSubsystem(const std::string& name = "CYLINDERSTRIP", const int layer = 0);

  //! destructor
  virtual ~PHG4CylinderStripSubsystem(void);

  //! init runwise stuff
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

// this method is used to check if it can be used as mothervolume
// Subsystems which can be mothervolume need to implement this 
// and return true
  virtual bool CanBeMotherSubsystem() const {return true;}

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4CylinderStripDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

};

#endif  // G4DETECTORS_PHG4CYLINDERSUBSYSTEM_H
