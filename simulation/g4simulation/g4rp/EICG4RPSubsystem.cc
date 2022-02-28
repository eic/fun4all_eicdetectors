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

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility> 

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
    std::string nodenameVirtSheet;

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

    nodenameVirtSheet = nodename + "_VirtSheet";

    PHG4HitContainer *rp_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
    PHG4HitContainer *rp_hitsVirtSheet = findNode::getClass<PHG4HitContainer>(topNode, nodenameVirtSheet);
    
    if( ! rp_hits )
    {
      dstNode->addNode(new PHIODataNode<PHObject>(rp_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
    }
    if( ! rp_hitsVirtSheet )
    {
	    dstNode->addNode(new PHIODataNode<PHObject>(rp_hits = new PHG4HitContainer(nodenameVirtSheet), nodenameVirtSheet, "PHObject"));
    }

    rp_hits->AddLayer(GetLayer());

    auto *tmp = new EICG4RPSteppingAction(this, m_Detector, GetParams());
    tmp->HitNodeName(nodename);
    tmp->HitNodeNameVirt(nodenameVirtSheet);
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

  std::string filename;
  std::string calibroot = std::string(getenv("CALIBRATIONROOT"));
  
  if ( calibroot.empty() )
  {    
    std::cout << "no CALIBRATIONROOT environment variable" << std::endl;
    gSystem->Exit(1);
  }

  filename = calibroot + "/RomanPots/RP_parameters_IP6.dat";

  set_default_string_param("parameter_file", filename);  
  set_default_double_param("FFenclosure_center", 2500);
  set_default_int_param("layerNumber", 1);
  set_default_double_param("detid", 0.);     //detector id
  set_default_int_param("lightyield", 0);
  set_default_int_param("use_g4steps", 0);
  set_default_double_param("tmin", NAN);
  set_default_double_param("tmax", NAN);

  set_default_double_param("Number_layers"   	 , 0);
  set_default_double_param("Sensor_size"   	 , 0);

  set_default_double_param("Layer1_pos_x"   	 , 0);
  set_default_double_param("Layer1_pos_y" 	 , 0);
  set_default_double_param("Layer1_pos_z" 	 , 0);
  set_default_double_param("Layer1_size_x" 	 , 0);
  set_default_double_param("Layer1_size_y" 	 , 0);
  set_default_double_param("Layer1_10sigma_x"    , 0);
  set_default_double_param("Layer1_10sigma_y"    , 0);
  set_default_double_param("Layer1_Si_size_z"    , 0);
  set_default_double_param("Layer1_Cu_size_z"    , 0);
  set_default_double_param("Layer1_rot_y"        , 0); 	

  set_default_double_param("Layer2_pos_x"   	 , 0);
  set_default_double_param("Layer2_pos_y" 	 , 0);
  set_default_double_param("Layer2_pos_z" 	 , 0);
  set_default_double_param("Layer2_size_x" 	 , 0);
  set_default_double_param("Layer2_size_y" 	 , 0);
  set_default_double_param("Layer2_10sigma_x"    , 0);
  set_default_double_param("Layer2_10sigma_y"    , 0);
  set_default_double_param("Layer2_Si_size_z"    , 0);
  set_default_double_param("Layer2_Cu_size_z"    , 0);
  set_default_double_param("Layer2_rot_y"        , 0);

}

//_______________________________________________________________________
void EICG4RPSubsystem::SetParameterFile(std::string &filename)
{
  set_string_param("parameter_file", filename );

  std::ifstream infile;
  std::string line;

  std::string paramFile = filename;   
  infile.open( paramFile );

  if(!infile.is_open()) {
    std::cout << "ERROR in EICG4RPDetector: Failed to open parameter file " << paramFile << std::endl;
   gSystem->Exit(1);
  }

  while( std::getline(infile, line) ) {
  
    std::string name;
    double value;
  
    std::istringstream iss( line );
  
    // skip comment lines
    if( line.find("#") != std::string::npos ) { continue; }
  
    if( !(iss >> name >> value) ) {
      std::cout << "Could not decode " << line << std::endl;
      gSystem->Exit(1);
    }

    set_double_param(name, value);

  }

  infile.close();

}
