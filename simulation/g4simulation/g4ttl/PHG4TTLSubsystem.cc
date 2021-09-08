#include "PHG4TTLSubsystem.h"
#include "PHG4TTLDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4Subsystem.h>       // for PHG4Subsystem
#include "PHG4TTLDisplayAction.h"
#include "PHG4TTLSteppingAction.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <sstream>

#include <set>  // for set
#include <sstream>
class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4TTLSubsystem::PHG4TTLSubsystem(const std::string& name)
  : PHG4DetectorSubsystem(name)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
  , superdetector("NONE")
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4TTLSubsystem::~PHG4TTLSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4TTLSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4TTLDisplayAction(Name(), showdetailed);
  // create detector
  m_Detector = new PHG4TTLDetector(this, topNode, GetParams(), Name());
  // m_Detector->geom = geom;
  m_Detector->SuperDetector(superdetector);
  m_Detector->OverlapCheck(CheckOverlap());

  // if (geom.GetNumActiveLayers())
  // {
  std::ostringstream nodename;
  if (superdetector != "NONE")
  {
    nodename << "G4HIT_" << superdetector;
  }
  else
  {
    nodename << "G4HIT_" << Name();
  }
  // create hit list
  PHG4HitContainer* block_hits = findNode::getClass<PHG4HitContainer>(
      topNode, nodename.str().c_str());
  if (!block_hits)
  {
    dstNode->addNode(new PHIODataNode<PHObject>(new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
  }
  // create stepping action
  m_SteppingAction = new PHG4TTLSteppingAction(m_Detector);
  m_Detector->SetSteppingAction(dynamic_cast<PHG4TTLSteppingAction*>(m_SteppingAction));
  // }
  return 0;
}

//_______________________________________________________________________
int PHG4TTLSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

//_______________________________________________________________________
PHG4Detector*
PHG4TTLSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4TTLSubsystem::SetDefaultParameters()
{
  set_default_int_param("isForward", 1);
  set_default_double_param("length", 100. * cm);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 375. * cm);
  set_default_double_param("offset_x", 150. * cm);
  set_default_double_param("rMin", 20. * cm);
  set_default_double_param("rMax", 220. * cm);
  set_default_double_param("etaMin", 1.);
  set_default_double_param("etaMax", 4.);
  set_default_double_param("polar_angle", 0.);
  set_default_double_param("total_thickness", 1. * cm);
  set_default_double_param("tSilicon", 1. * um);
  set_default_double_param("resLGAD", 500e-4 / sqrt(12));
  // set_default_double_param("wls_dw", 0.3);
  // set_default_double_param("support_dw", 0.2);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  // set_default_double_param("thickness_absorber", 2.);
  // set_default_double_param("thickness_scintillator", 0.231);
  // set_default_string_param("scintillator", "G4_POLYSTYRENE");
  // set_default_string_param("absorber", "G4_Fe");
  // set_default_string_param("support", "G4_Fe");
  return;
}
