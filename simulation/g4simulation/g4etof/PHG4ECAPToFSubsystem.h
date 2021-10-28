#ifndef G4DETECTORS_PHG4ECAPToFSUBSYSTEM_H
#define G4DETECTORS_PHG4ECAPToFSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <array>   // for array
#include <string>  // for string

class PHCompositeNode;
class PHG4ECAPToFDetector;
class PHG4Detector;
class PHG4SteppingAction;

class PHG4ECAPToFSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4ECAPToFSubsystem(const std::string& name = "ECAPToF", const int layer = 0);

  //! destructor
  ~PHG4ECAPToFSubsystem(void) override;

  int InitRunSubsystem(PHCompositeNode*) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
*/
  int process_event(PHCompositeNode*) override;

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const override;
  PHG4SteppingAction* GetSteppingAction(void) const override { return m_SteppingAction; }

  // this method is used to check if it can be used as mothervolume
  // Subsystems which can be mothervolume need to implement this
  // and return true
  //bool CanBeMotherSubsystem() const override { return true; }

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4ECAPToFDetector* m_Detector /*= nullptr*/;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction /*= nullptr*/;

  //bool m_SaveAllHitsFlag = false;
};
#endif  // G4DETECTORS_PHG4ECAPToFSUBSYSTEM_H
