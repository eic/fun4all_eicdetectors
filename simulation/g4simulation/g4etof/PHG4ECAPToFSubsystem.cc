#include "PHG4ECAPToFSubsystem.h"
#include "PHG4ECAPToFDetector.h"
#include "PHG4ECAPToFSteppingAction.h"

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
#include <iostream>  // for operator<<, basic_ostream, std::endl
#include <sstream>

class PHG4Detector;
class PHG4SteppingAction;

PHG4ECAPToFSubsystem::PHG4ECAPToFSubsystem(const std::string &na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
{
  InitializeParameters();
}

PHG4ECAPToFSubsystem::~PHG4ECAPToFSubsystem()
{
}

//_______________________________________________________________________
int PHG4ECAPToFSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  // create detector
  m_Detector = new PHG4ECAPToFDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  /*
  string nodename;
  nodename = "G4HIT_" + SuperDetector();
  string nodes;
  
  if (GetParams()->get_int_param("active"))
  {
  PHG4HitContainer *ECAPToF_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
  
  if (!ECAPToF_hits)
  {
  //dstNode->addNode(new PHIODataNode<PHObject>(ECAPToF_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
  ECAPToF_hits = new PHG4HitContainer(nodename);
  PHIODataNode<PHObject>* hitNode = new PHIODataNode<PHObject>(ECAPToF_hits, nodename, "PHObject");
  dstNode->addNode(hitNode);
  // nodes.insert(nodename);
  }
  */

  std::set<std::string> nodes;
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }

    std::ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    nodes.insert(nodename.str());
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodename.str("");
      if (SuperDetector() != "NONE")
      {
        nodename << "G4HIT_ACTIVEGAS_" << SuperDetector();
      }
      else
      {
        nodename << "G4HIT_ACTIVEGAS_" << Name();
      }
      nodes.insert(nodename.str());
    }

    for (std::string node : nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node.c_str());
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node.c_str(), "PHObject"));
      }
    }
    //Stepping action
    /*
	auto *tmp = new PHG4ECAPToFSteppingAction(this, m_Detector, GetParams());
	tmp->HitNodeName(nodename);
	m_SteppingAction = tmp;
	}
	else if (GetParams()->get_int_param("blackhole"))
  {
  m_SteppingAction = new PHG4ECAPToFSteppingAction(this, m_Detector, GetParams());
  }
  
  if (m_SteppingAction)
  {
  (dynamic_cast<PHG4ECAPToFSteppingAction *>(m_SteppingAction))->SaveAllHits(m_SaveAllHitsFlag);
  }
  */
    m_SteppingAction = new PHG4ECAPToFSteppingAction(m_Detector, GetParams());
  }
  else
  {
    // if this is a black hole it does not have to be active
    if (GetParams()->get_int_param("blackhole"))
    {
      m_SteppingAction = new PHG4ECAPToFSteppingAction(m_Detector, GetParams());
    }
  }

  return 0;
}

int PHG4ECAPToFSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4ECAPToFSubsystem::SetDefaultParameters()
{
  set_default_int_param("n_fgas_layer", 6);
  set_default_int_param("n_bgas_layer", 6);
  set_default_double_param("gas_gap", 0.04);
  set_default_double_param("glass_thick", 0.022);
  set_default_double_param("Carbon_thick", 0.01);
  set_default_double_param("pcb_thick", 0.04);
  set_default_double_param("cu_thick", 0.01);
  set_default_double_param("honeycomb_thick", 0.075);
  set_default_double_param("mylar_thick", 0.01);
  set_default_double_param("Rin", 5.);
  set_default_double_param("Rout", 168.0);
  set_default_double_param("z_begin", 72.0);
  set_default_int_param("use_g4steps", 1);

  // place holder, will be replaced by world material if not set by other means (macro)
  set_default_string_param("material", "G4_AIR");
}

PHG4Detector *PHG4ECAPToFSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4ECAPToFSubsystem::Print(const std::string &what) const
{
  std::cout << Name() << " Parameters: " << std::endl;
  if (!BeginRunExecuted())
  {
    std::cout << "Need to execute BeginRun() before parameter printout is meaningful" << std::endl;
    std::cout << "To do so either run one or more events or on the command line execute: " << std::endl;
    std::cout << "Fun4AllServer *se = Fun4AllServer::instance();" << std::endl;
    std::cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << std::endl;
    std::cout << "g4->InitRun(se->topNode());" << std::endl;
    std::cout << "PHG4ECAPToFSubsystem *trd = (PHG4ECAPToFSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << std::endl;
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
