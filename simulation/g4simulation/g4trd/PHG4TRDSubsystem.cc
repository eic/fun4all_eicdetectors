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

#include <boost/foreach.hpp>

#include <cmath>     // for NAN
#include <iostream>  // for operator<<, basic_ostream, endl
#include <sstream>

class PHG4Detector;
class PHG4SteppingAction;

using namespace std;

PHG4TRDSubsystem::PHG4TRDSubsystem(const std::string &na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr),
    m_Detector(nullptr),
    m_SteppingAction(nullptr)
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
/*
  string nodename;
  nodename = "G4HIT_" + SuperDetector();
  string nodes;
  
  if (GetParams()->get_int_param("active"))
  {
  PHG4HitContainer *TRD_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
  
  if (!TRD_hits)
  {
  //dstNode->addNode(new PHIODataNode<PHObject>(TRD_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
  TRD_hits = new PHG4HitContainer(nodename);
  PHIODataNode<PHObject>* hitNode = new PHIODataNode<PHObject>(TRD_hits, nodename, "PHObject");
  dstNode->addNode(hitNode);
  // nodes.insert(nodename);
  }
*/
  
  set<string> nodes;
  if (GetParams()->get_int_param("active"))
    {
      PHNodeIterator dstIter(dstNode);
      PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
	{
	  DetNode = new PHCompositeNode(SuperDetector());
	  dstNode->addNode(DetNode);
	}
      
      ostringstream nodename;
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
      
      BOOST_FOREACH (string node, nodes)
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
	auto *tmp = new PHG4TRDSteppingAction(this, m_Detector, GetParams());
	tmp->HitNodeName(nodename);
	m_SteppingAction = tmp;
	}
	else if (GetParams()->get_int_param("blackhole"))
  {
  m_SteppingAction = new PHG4TRDSteppingAction(this, m_Detector, GetParams());
  }
  
  if (m_SteppingAction)
  {
  (dynamic_cast<PHG4TRDSteppingAction *>(m_SteppingAction))->SaveAllHits(m_SaveAllHitsFlag);
  }
*/
      m_SteppingAction = new PHG4TRDSteppingAction(m_Detector, GetParams());
    }
  else
    {
// if this is a black hole it does not have to be active
      if (GetParams()->get_int_param("blackhole"))
	{
	  m_SteppingAction = new PHG4TRDSteppingAction(m_Detector, GetParams());
	}
    }
 
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
 
 PHG4Detector *  PHG4TRDSubsystem::GetDetector(void) const
 {
   return m_Detector;
 }
 
     
void PHG4TRDSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  if (!BeginRunExecuted())
    {
      cout << "Need to execute BeginRun() before parameter printout is meaningful" << endl;
      cout << "To do so either run one or more events or on the command line execute: " << endl;
      cout << "Fun4AllServer *se = Fun4AllServer::instance();" << endl;
      cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << endl;
      cout << "g4->InitRun(se->topNode());" << endl;
      cout << "PHG4TRDSubsystem *trd = (PHG4TRDSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << endl;
      cout << "trd->Print()" << endl;
      return;
    }
  GetParams()->Print();
  if (m_SteppingAction)
    {
      m_SteppingAction->Print(what);
    }
  return;
}
     
