//____________________________________________________________________________..
//
// This is the interface to the framework. You only need to define the parameters
// you use for your detector in the SetDefaultParameters() method here
// The place to do this is marked by //implement your own here//
// The parameters have no units, they need to be converted in the
// EICG4RPDetector::ConstructMe() method
// but the convention is as mentioned cm and deg
//____________________________________________________________________________..
//
#include "EICG4RPSubsystem.h"

#include "EICG4RPDetector.h"
#include "EICG4RPSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

//_______________________________________________________________________
EICG4RPSubsystem::EICG4RPSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
{
  // call base class method which will set up parameter infrastructure
  // and call our SetDefaultParameters() method
  InitializeParameters();
}
//_______________________________________________________________________
int EICG4RPSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  // create detector
  m_Detector = new EICG4RPDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    std::string nodename;

    if (SuperDetector() != "NONE")
    {
      // create super detector subnodes
      PHNodeIterator iter_dst(dstNode);
      PHCompositeNode *superSubNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(superSubNode);
      }
      dstNode = superSubNode;

      nodename = "G4HIT_" + SuperDetector();
    }

    else
    {
      nodename = "G4HIT_" + Name();
    }
    PHG4HitContainer *rp_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
    if (!rp_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(rp_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
    }
    rp_hits->AddLayer(GetLayer());
    auto *tmp = new EICG4RPSteppingAction(this, m_Detector, GetParams());
    tmp->HitNodeName(nodename);
    m_SteppingAction = tmp;
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new EICG4RPSteppingAction(this, m_Detector, GetParams());
  }
  if (m_SteppingAction)
  {
    (dynamic_cast<EICG4RPSteppingAction *>(m_SteppingAction))->SaveAllHits(m_SaveAllHitsFlag);
  }
  return 0;
}
//_______________________________________________________________________
int EICG4RPSubsystem::process_event(PHCompositeNode *topNode)
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
void EICG4RPSubsystem::Print(const std::string &what) const
{
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *EICG4RPSubsystem::GetDetector(void) const
{
  return m_Detector;
}

//_______________________________________________________________________
void EICG4RPSubsystem::SetDefaultParameters()
{
  // sizes are in cm
  // angles are in deg
  // units should be converted to G4 units when used
  //implement your own here//
  set_default_double_param("place_x", 0.);   //subdetector position
  set_default_double_param("place_y", 0.);   //subdetector position
  set_default_double_param("place_z", 0.);   //subdetector position
  set_default_double_param("hole_x", 10.0);  //beam pipe cut off in the detector volume
  set_default_double_param("hole_y", 5.0);   //beam pipe cut off in the detector volume
  set_default_double_param("pipe_x", 0);     //beam pipe position
  set_default_double_param("pipe_y", 0.);    //beam pipe position
  set_default_double_param("pipe_z", 0.);    //beam pipe position
  set_default_double_param("rot_y", 0.025);  //subdetector rotation
  set_default_double_param("rp_x", 25.);     //detector outer radiues
  set_default_double_param("rp_y", 10.);     //detector outer radiues
  set_default_double_param("length", .1);    //detector length
  set_default_double_param("detid", 0.);     //detector id
  set_default_int_param("lightyield", 0);
  set_default_int_param("use_g4steps", 0);
  set_default_double_param("tmin", NAN);
  set_default_double_param("tmax", NAN);
  set_default_string_param("material", "G4_PbWO4");  //detector material
}