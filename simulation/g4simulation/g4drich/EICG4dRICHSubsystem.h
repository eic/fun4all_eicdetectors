#ifndef DRICHSUBSYSTEM_H
#define DRICHSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

class PHCompositeNode;
class PHG4Detector;
class EICG4dRICHDetector;
class PHG4SteppingAction;

/**
 * \brief Detector Subsystem module
 *
 * The detector is constructed and registered via EICG4dRICHDetector
 *
 *
 * \see EICG4dRICHDetector
 * \see EICG4dRICHSubsystem
 *
 */
class EICG4dRICHSubsystem : public PHG4DetectorSubsystem 
{
  public:
    //! constructor
    EICG4dRICHSubsystem(const std::string &name = "EICG4dRICH");

    //! destructor
    virtual ~EICG4dRICHSubsystem() {}

    /*!
    creates relevant hit nodes that will be populated by the stepping action and
    stored in the output DST
    */
    int InitRunSubsystem(PHCompositeNode *) override;

    //! event processing
    /*!
    get all relevant nodes from top nodes (namely hit list)
    and pass that to the stepping action
    */
    int process_event(PHCompositeNode *) override;

    //! accessors (reimplemented)
    PHG4Detector *GetDetector() const override;

    PHG4SteppingAction *GetSteppingAction() const override { return m_SteppingAction; }

    void SetGeometryFile(const std::string &fileName) { m_geoFile = fileName; }   
 
    //! Print info (from SubsysReco)
    void Print(const std::string &what = "ALL") const override;

  protected:
    // \brief Set default parameter values
    void SetDefaultParameters() override;

  private:
    //! detector construction
    /*! derives from PHG4Detector */
    EICG4dRICHDetector *m_Detector;

    //! particle tracking "stepping" action
    /*! derives from PHG4SteppingActions */
    PHG4SteppingAction *m_SteppingAction;
   
    std::string m_geoFile;
};

#endif // DRICHSUBSYSTEM_H
