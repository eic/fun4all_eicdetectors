#include "PHG4TRDSubsystem.h"
#include "PHG4TRDDetector.h"
#include "PHG4TRDSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <cmath>     // for NAN
#include <iostream>  // for operator<<, basic_ostream, endl
#include <sstream>

class PHG4Detector;
class PHG4SteppingAction;

PHG4TRDSubsystem::PHG4TRDSubsystem(const std::string &na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
{
  InitializeParameters();
}

PHG4TRDSubsystem::~PHG4TRDSubsystem()
{
}

//_______________________________________________________________________
int PHG4TRDSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  // create detector
  m_Detector = new PHG4TRDDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  std::string nodename;
  nodename = "G4HIT_" + SuperDetector();
  //std::string nodes;

  if (GetParams()->get_int_param("active"))
  {
    PHG4HitContainer *TRD_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);

    if (!TRD_hits)
    {
      //dstNode->addNode(new PHIODataNode<PHObject>(TRD_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
      TRD_hits = new PHG4HitContainer(nodename);
      PHIODataNode<PHObject> *hitNode = new PHIODataNode<PHObject>(TRD_hits, nodename, "PHObject");
      dstNode->addNode(hitNode);
      // nodes.insert(nodename);
    }

    //Stepping action
    auto *tmp = new PHG4TRDSteppingAction(this, m_Detector, GetParams());
    tmp->HitNodeName(nodename);
    m_SteppingAction = tmp;
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new PHG4TRDSteppingAction(this, m_Detector, GetParams());
  }
  /*
  if (m_SteppingAction)
    {
      (dynamic_cast<PHG4TRDSteppingAction *>(m_SteppingAction))->SaveAllHits(m_SaveAllHitsFlag);
    }
    */
  return 0;
}

int PHG4TRDSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4TRDSubsystem::SetDefaultParameters()
{
  set_default_double_param("ThicknessZ", 35.);
  set_default_double_param("RIn", 20.);
  set_default_double_param("ROut", 200.);
  set_default_double_param("PosZ", 30.);
  set_default_double_param("det_RIn", 20.);
  set_default_double_param("det_ROut", 200.);
  set_default_int_param("use_g4steps", 1);

  // place holder, will be replaced by world material if not set by other means (macro)
  set_default_string_param("material", "G4_AIR");
}

PHG4Detector *PHG4TRDSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4TRDSubsystem::Print(const std::string &what) const
{
  std::cout << Name() << " Parameters: " << std::endl;
  if (!BeginRunExecuted())
  {
    std::cout << "Need to execute BeginRun() before parameter printout is meaningful" << std::endl;
    std::cout << "To do so either run one or more events or on the command line execute: " << std::endl;
    std::cout << "Fun4AllServer *se = Fun4AllServer::instance();" << std::endl;
    std::cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << std::endl;
    std::cout << "g4->InitRun(se->topNode());" << std::endl;
    std::cout << "PHG4TRDSubsystem *trd = (PHG4TRDSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << std::endl;
    std::cout << "trd->Print()" << std::endl;
    return;
  }
  GetParams()->Print();
  if (m_SteppingAction)
  {
    m_SteppingAction->Print(what);
  }
  return;
}
