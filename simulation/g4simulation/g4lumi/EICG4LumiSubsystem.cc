//____________________________________________________________________________..
//
// This is the interface to the framework. You only need to define the parameters
// you use for your detector in the SetDefaultParameters() method here
// The place to do this is marked by //implement your own here//
// The parameters have no units, they need to be converted in the
// EICG4LumiDetector::ConstructMe() method
// but the convention is as mentioned cm and deg
//____________________________________________________________________________..
//
#include "EICG4LumiSubsystem.h"

#include "EICG4LumiDetector.h"
#include "EICG4LumiSteppingAction.h"

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
EICG4LumiSubsystem::EICG4LumiSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
{
  // call base class method which will set up parameter infrastructure
  // and call our SetDefaultParameters() method
  InitializeParameters();
}

//_______________________________________________________________________
int EICG4LumiSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  // create detector
  m_Detector = new EICG4LumiDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  
     std::string nodename_prefix;
     std::string nodename_photonCAL;
     std::string nodename_virtExitWindow;
   
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
   
       nodename_prefix = "G4HIT_" + SuperDetector();
    }
 
    else
    {
       nodename_prefix = "G4HIT_" + Name();
    }
    
    nodename_photonCAL = nodename_prefix + "_photonCAL";
    nodename_virtExitWindow = nodename_prefix + "_virtExitWindow";
    
    // Grab hit container if it exists
    PHG4HitContainer *PhotonCAL_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename_photonCAL);
    PHG4HitContainer *VirtExitWindow_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename_virtExitWindow);
    // Create it if it doesn't exist
    if( ! PhotonCAL_hits )
    {
      dstNode->addNode( new PHIODataNode<PHObject>( PhotonCAL_hits = new PHG4HitContainer(nodename_photonCAL), nodename_photonCAL, "PHObject"));
    }
    if( ! VirtExitWindow_hits )
    {
      dstNode->addNode( new PHIODataNode<PHObject>( VirtExitWindow_hits = new PHG4HitContainer(nodename_virtExitWindow), nodename_virtExitWindow, "PHObject"));
    }

    PhotonCAL_hits->AddLayer( GetLayer() );
    VirtExitWindow_hits->AddLayer( GetLayer() );
    
    // Create stepping action
    auto *tmp = new EICG4LumiSteppingAction( this, m_Detector, GetParams() );
    tmp->HitNodeNameVirt( nodename_virtExitWindow) ;
    tmp->HitNodeName( nodename_photonCAL );
    
    m_SteppingAction = tmp;

  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new EICG4LumiSteppingAction(this, m_Detector, GetParams());
  }

  if (m_SteppingAction)
  {
    (dynamic_cast<EICG4LumiSteppingAction *>(m_SteppingAction))->SaveAllHits(m_SaveAllHitsFlag);
  }

    return 0;
}
//_______________________________________________________________________
int EICG4LumiSubsystem::process_event(PHCompositeNode *topNode)
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
void EICG4LumiSubsystem::Print(const string &what) const
{
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *EICG4LumiSubsystem::GetDetector(void) const
{
  return m_Detector;
}

//_______________________________________________________________________
void EICG4LumiSubsystem::SetDefaultParameters()
{
 
// defaults are in cm and rad
  	set_default_double_param("FBenclosure_center", 0.);
  	set_default_string_param("parameter_file", "");
  	set_default_int_param("layerNumber", 1);
  	set_default_double_param("detid", 0.);     //detector id
  	set_default_int_param("lightyield", 0);
  	set_default_int_param("use_g4steps", 0);
  	set_default_double_param("tmin", NAN);
  	set_default_double_param("tmax", NAN);
  //
  	set_default_double_param("LumiWin_X", 0.);
  	set_default_double_param("LumiWin_Y", 0.);
  	set_default_double_param("LumiWin_Z", -1850);
  	set_default_double_param("LumiWin_Tilt", 0.25);
  	set_default_double_param("LumiWin_Thickness", 0.26);
  	set_default_double_param("LumiWin_Height", 7.4);
  	set_default_double_param("LumiWin_Length", 29);
  	set_default_string_param("LumiWin_Material", "G4_Al");
  //
  	set_default_double_param("LumiMag_Z", -2820);
  	set_default_double_param("LumiMag_innerR", 10);
  	set_default_double_param("LumiMag_outerR", 16);
  	set_default_double_param("LumiMag_DZ", 60);
  	set_default_double_param("LumiMag_B", 0.37);
  	set_default_string_param("LumiMag_VesselMaterial", "G4_Fe");
//
  	set_default_double_param("LumiSpec_Z", -3640);
	set_default_double_param("LumiSpec_XY", 20);
	set_default_double_param("LumiSpec_DZ", 35);
	set_default_double_param("LumiSpec_YS", 6);
//
	set_default_double_param("LumiPhotonCAL_Z", -3700);
	set_default_double_param("LumiPhotonCAL_XY", 16);
	set_default_double_param("LumiPhotonCAL_DZ", 35);	

}
