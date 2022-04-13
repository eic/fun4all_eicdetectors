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

  // maximum of 4 layers for the IP8 config, only layers 1 and 2 used for IP6
  for( int l=1; l <= 4; l++ ) {
    set_default_double_param(Form("Layer%i_pos_x", l)   	, 0);
    set_default_double_param(Form("Layer%i_pos_y", l) 	        , 0);
    set_default_double_param(Form("Layer%i_pos_z", l) 	        , 0);
    set_default_double_param(Form("Layer%i_size_x", l) 	        , 0);
    set_default_double_param(Form("Layer%i_size_y", l) 	        , 0);
    set_default_double_param(Form("Layer%i_beamHoleHW_x", l)    , 0);
    set_default_double_param(Form("Layer%i_beamHoleHW_y", l)    , 0);
    set_default_double_param(Form("Layer%i_Si_size_z", l)       , 0);
    set_default_double_param(Form("Layer%i_Cu_size_z", l)       , 0);
    set_default_double_param(Form("Layer%i_rot_y", l)           , 0); 	
  }     
  
}

//_______________________________________________________________________
void EICG4RPSubsystem::SetParametersFromFile( std::string filename )
{
  set_string_param("parameter_file", filename );

  std::ifstream infile;
  std::string line;

  std::string paramFile = filename;   
  infile.open( paramFile );

  if(!infile.is_open()) {
    std::cout << "ERROR in EICG4RPSubsystem: Failed to open parameter file " << paramFile << std::endl;
    gSystem->Exit(1);
  }

  while( std::getline(infile, line) ) {

    std::string name;
    double value;

    std::istringstream iss( line );

    // skip comment lines
    if( line.find("#") != std::string::npos ) { continue; }

    // grab the line
    if( !(iss >> name >> value) ) {
      std::cout << "Could not decode " << line << std::endl;
      gSystem->Exit(1);
    }

    // simply set all parameters NOT related to the beam hole
    if( name.find("beamHoleHW") == std::string::npos ) {
      set_double_param( name, value );
    }
    else { // parse beam hole size based on chosen beam energy and config

      std::string prefix = name.substr( 0, name.find("__") ); // LayerN_beamHoleHW_x(or)y
      std::string expected_string = prefix;

      if( m_beamProfile.find("eA") != std::string::npos ) { // Heavy ions
        expected_string += "__eA";

        double ionE_diff_1 = fabs( m_ionE - 110 );
        double ionE_diff_2 = fabs( m_ionE - 41 );
        double elecE_diff_1 = fabs( m_elecE - 18 );
        double elecE_diff_2 = fabs( m_elecE - 10 );
        double elecE_diff_3 = fabs( m_elecE - 5 );

        if( ionE_diff_2 < ionE_diff_1 ) {
          expected_string += "_41x5";
        }
        else {
          if( elecE_diff_3 < elecE_diff_2 )      { expected_string += "_110x5"; }
          else if( elecE_diff_2 < elecE_diff_1 ) { expected_string += "_110x10"; }
          else                                   { expected_string += "_110x18"; }
        }
      }
      else { // protons
        if( m_beamProfile.find("high-acceptance") != std::string::npos ) { 
          expected_string += "__epHA"; 
        }
        else { 
          expected_string += "__epHD"; 
        }

        double ionE_diff_1 = fabs( m_ionE - 275 );
        double ionE_diff_2 = fabs( m_ionE - 100 );
        double ionE_diff_3 = fabs( m_ionE - 41 );
        double elecE_diff_1 = fabs( m_elecE - 18 );
        double elecE_diff_2 = fabs( m_elecE - 10 );
        double elecE_diff_3 = fabs( m_elecE - 5 );

        if( ionE_diff_3 < ionE_diff_2 ) { // 41 GeV protons
          expected_string += "_41x5";
        }
        else if( ionE_diff_2 < ionE_diff_1 ) { // 100 GeV protons
          if( elecE_diff_3 < elecE_diff_2 )  { expected_string += "_100x5"; }
          else                               { expected_string += "_100x10"; }
        }
        else {
          if( elecE_diff_2 < elecE_diff_1 )  { expected_string += "_275x10"; }
          else                               { expected_string += "_275x18"; }
        }
      }

      // skip undesired beam hole configuration
      if( name.compare( expected_string ) != 0 ) { continue; }
      
      set_double_param( prefix, value);
    }

  }

  infile.close();

}
