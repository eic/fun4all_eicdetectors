// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICG4LUMISUBSYSTEM_H
#define EICG4LUMISUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

class PHCompositeNode;
class PHG4Detector;
class EICG4LumiDetector;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   *
   * The detector is constructed and registered via EICG4LumiDetector
   *
   *
   * \see EICG4LumiDetector
   * \see EICG4LumiSubsystem
   *
   */
class EICG4LumiSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  EICG4LumiSubsystem(const std::string& name = "EICG4Lumi", const int layer = 0);

  //! destructor
  virtual ~EICG4LumiSubsystem() {}

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

  PHG4SteppingAction* GetSteppingAction() const override { return m_SteppingAction; }  //StepAc

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

  bool CanBeMotherSubsystem() const override {return true;}

  void SaveAllHits(bool i = true){ m_SaveAllHitsFlag = i;}  
  //void SetTowerMappingFile(const std::string &filename);

  void SetParameterFile( std::string &filename ) { set_string_param("parameter_file", filename ); } //////

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  EICG4LumiDetector  *m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  
  PHG4SteppingAction *m_SteppingAction; //StepAc

  bool m_SaveAllHitsFlag = false;
  //std::string mappingfile_;

};

#endif // EICG4LUMISUBSYSTEM_H



