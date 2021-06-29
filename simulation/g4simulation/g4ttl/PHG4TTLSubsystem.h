// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TTLSUBSYSTEM_H
#define G4DETECTORS_PHG4TTLSUBSYSTEM_H

// #include "PHG4TTLDetector.h"

#include <g4detectors/PHG4DetectorSubsystem.h>  // for PHG4DetectorSubsystem

#include <string>                   // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4TTLDetector;
class PHG4SteppingAction;

class PHG4TTLSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4TTLSubsystem(const std::string& name = "Sector");

  //! destructor
  ~PHG4TTLSubsystem() override;

  //! init
  /*!
   creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
   reates the stepping action and place it on the node tree, under "ACTIONS" node
   creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
   */
  int InitRunSubsystem(PHCompositeNode *) override;

  //! event processing
  /*!
   get all relevant nodes from top nodes (namely hit list)
   and pass that to the stepping action
   */
  int process_event(PHCompositeNode*) override;

  //! accessors (reimplemented)
  PHG4Detector*
  GetDetector(void) const override;
  PHG4SteppingAction* GetSteppingAction(void) const override { return m_SteppingAction; }

  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

  void
  SuperDetector(const std::string& name)
  {
    superdetector = name;
  }

  // //! geometry manager PHG4TTL::Sector_Geometry
  // PHG4TTL::Sector_Geometry&
  // get_geometry()
  // {
  //   return geom;
  // }

  // //! geometry manager PHG4TTL::Sector_Geometry
  // void
  // set_geometry(const PHG4TTL::Sector_Geometry& g)
  // {
  //   geom = g;
  // }

  /** Set level of detail for display
   */
  void SetDetailed(bool b){showdetailed = b;}

  
 private:
  void SetDefaultParameters() override;
  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4TTLDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  std::string superdetector;
  bool showdetailed = false;

  // PHG4TTL::Sector_Geometry geom;
};

#endif
