// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4BwdSUBSYSTEM_H
#define EICG4BwdSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

class PHCompositeNode;
class PHG4Detector;
class EICG4BwdDetector;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   *
   * The detector is constructed and registered via EICG4B0Detector
   *
   *
   * \see EICG4BwdDetector
   * \see EICG4BwdSubsystem
   *
   */
class EICG4BwdSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  EICG4BwdSubsystem(const std::string& name = "EICG4Bwd", const int layer = 0);

  //! destructor
  virtual ~EICG4BwdSubsystem() {}

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

  bool CanBeMotherSubsystem() const override {return true;}

  void SaveAllHits(bool i = true){ m_SaveAllHitsFlag = i;}

  void SetTowerMappingFile(const std::string &filename);

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  EICG4BwdDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;
  
  bool m_SaveAllHitsFlag = false;

  std::string mappingfile_;
};

#endif  // EICG4BwdSUBSYSTEM_H
