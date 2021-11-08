//____________________________________________________________________________..
//
// This is the interface to the framework. You only need to define the parameters
// you use for your detector in the SetDefaultParameters() method here
// The place to do this is marked by //implement your own here//
// The parameters have no units, they need to be converted in the
// EICG4B0Detector::ConstructMe() method
// but the convention is as mentioned cm and deg
//____________________________________________________________________________..
//
#include "EICG4B0Subsystem.h"

#include "EICG4B0Detector.h"
#include "EICG4B0SteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

using namespace std;

//_______________________________________________________________________
EICG4B0Subsystem::EICG4B0Subsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , mappingfile_("")
{
  // call base class method which will set up parameter infrastructure
  // and call our SetDefaultParameters() method
  InitializeParameters();
}
//_______________________________________________________________________
int EICG4B0Subsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  // create detector
  m_Detector = new EICG4B0Detector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  
    string nodename;
    string geonode;

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
    PHG4HitContainer *b0_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
    if (!b0_hits)
    {
       dstNode->addNode(new PHIODataNode<PHObject>(b0_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
    }
    b0_hits->AddLayer(GetLayer());
    auto *tmp = new EICG4B0SteppingAction(this, m_Detector, GetParams());
    tmp->HitNodeName(nodename);
    m_SteppingAction = tmp;
   }
   else if (GetParams()->get_int_param("blackhole"))
   {
     m_SteppingAction = new EICG4B0SteppingAction(this, m_Detector, GetParams());
   }
   if (m_SteppingAction)
   {
      (dynamic_cast<EICG4B0SteppingAction *>(m_SteppingAction))->SaveAllHits(m_SaveAllHitsFlag);
   }
  return 0;
}
//_______________________________________________________________________
int EICG4B0Subsystem::process_event(PHCompositeNode *topNode)
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
void EICG4B0Subsystem::Print(const string &what) const
{
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *EICG4B0Subsystem::GetDetector(void) const
{
  return m_Detector;
}

//_______________________________________________________________________
void EICG4B0Subsystem::SetDefaultParameters()
{
  // sizes are in cm
  // angles are in deg
  // units should be converted to G4 units when used
  //implement your own here//
  set_default_double_param("place_x", 0.);          //subdetector position
  set_default_double_param("place_y", 0.);          //subdetector position
  set_default_double_param("place_z", 0.);          //subdetector position
  set_default_double_param("pipe_ir", 2.8);         //beam pipe inner radius (for future implementation)
  set_default_double_param("pipe_or", 3.05);        //beam pipe outer raidus (for future implementation)
  set_default_double_param("pipe_hole", 5.0);       //beam pipe cut off radius in the detector volume
  set_default_double_param("pipe_x", -3.4);         //beam pipe position
  set_default_double_param("pipe_y", 0.);           //beam pipe position
  set_default_double_param("pipe_z", 0.);           //beam pipe position
  set_default_double_param("rot_y", 0.);            //subdetector rotation
  set_default_double_param("outer_radius", 2.);     //detector outer radiues
  set_default_double_param("d_radius", 5.);         //packman cutoff size
  set_default_double_param("length", 10.);          //detector length
  set_default_double_param("startAngle", 0.);       //start Angle for packman cutoff
  set_default_double_param("spanningAngle", 360.);  //spanning Angle of the detector (for packman cutoff)
  set_default_double_param("detid", 0.);            //detector id
  set_default_int_param("ispipe", 0);               //pipe or detector (for future implementation)
  set_default_int_param("lightyield", 0);
  set_default_int_param("use_g4steps", 0);
  set_default_double_param("tower_size", 2.);
  set_default_double_param("readout_size", 2.);
  set_default_double_param("tmin", NAN);
  set_default_double_param("tmax", NAN);
  set_default_string_param("material", "G4_PbWO4");  //detector material
}

void EICG4B0Subsystem::SetTowerMappingFile(const std::string& filename)
{
  mappingfile_ = filename;
}
