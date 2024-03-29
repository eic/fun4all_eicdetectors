#include "PHG4HybridHomogeneousCalorimeterSubsystem.h"
#include "PHG4HybridHomogeneousCalorimeterDetector.h"
#include "PHG4HybridHomogeneousCalorimeterDisplayAction.h"
#include "PHG4HybridHomogeneousCalorimeterSteppingAction.h"
#include "PHG4ProjCrystalCalorimeterDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4Subsystem.h>       // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <iostream>  // for operator<<, ostrin...
#include <sstream>

class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4HybridHomogeneousCalorimeterSubsystem::PHG4HybridHomogeneousCalorimeterSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4HybridHomogeneousCalorimeterSubsystem::~PHG4HybridHomogeneousCalorimeterSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4HybridHomogeneousCalorimeterSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  // create display settings before detector
  m_DisplayAction = new PHG4HybridHomogeneousCalorimeterDisplayAction(Name());
  // create detector
  if (Verbosity() > 1)
  {
    cout << "PHG4HybridHomogeneousCalorimeterSubsystem::InitRun - use PHG4HybridHomogeneousCalorimeterDetector" << endl;
  }
  m_Detector = new PHG4HybridHomogeneousCalorimeterDetector(this, topNode, GetParams(), Name());

  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->DetectorId(GetLayer());
  m_Detector->DoFullLightProp(_do_lightpropagation);

  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

    // create hit output node
    string nodename = "G4HIT_";
    if (SuperDetector() != "NONE")
    {
      PHNodeIterator iter_dst(dstNode);
      PHCompositeNode* superSubNode = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(superSubNode);
      }
      dstNode = superSubNode;
      nodename += SuperDetector();
    }
    else
    {
      nodename += Name();
    }

    PHG4HitContainer* crystal_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
    if (!crystal_hits)
    {
      crystal_hits = new PHG4HitContainer(nodename);
      PHIODataNode<PHObject>* hitNode = new PHIODataNode<PHObject>(crystal_hits, nodename, "PHObject");
      dstNode->addNode(hitNode);
    }

    if (GetParams()->get_int_param("absorberactive"))
    {
      string absnodename = "G4HIT_ABSORBER_";
      if (SuperDetector() != "NONE")
      {
        absnodename += SuperDetector();
      }
      else
      {
        absnodename += Name();
      }

      PHG4HitContainer* absorber_hits = findNode::getClass<PHG4HitContainer>(topNode, absnodename);
      if (!absorber_hits)
      {
        absorber_hits = new PHG4HitContainer(absnodename);
        PHIODataNode<PHObject>* abshitNode = new PHIODataNode<PHObject>(absorber_hits, absnodename, "PHObject");
        dstNode->addNode(abshitNode);
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4HybridHomogeneousCalorimeterSteppingAction(m_Detector, GetParams());
  }
  return 0;
}

//_______________________________________________________________________
int PHG4HybridHomogeneousCalorimeterSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4HybridHomogeneousCalorimeterSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4HybridHomogeneousCalorimeterSubsystem::SetDefaultParameters()
{
  // values in cm and degrees
  set_default_int_param("projective", 0);
  set_default_int_param("carbon_frame_style", 0);
  set_default_double_param("carbon_frame_depth", 1.0);
  set_default_double_param("carbon_face_lip", 0.1);
  set_default_double_param("reflective_foil_thickness", 0.0);
  set_default_double_param("tedlar_thickness", 0.0);
  set_default_double_param("sensor_thickness", 0.0);
  set_default_double_param("sensor_dimension", 0.0);
  set_default_int_param("sensor_count", 0);

  set_default_double_param("crystal_dx", 2.);
  set_default_double_param("crystal_dy", 2.);
  set_default_double_param("crystal_dz", 18.);
  set_default_double_param("dz", 18.);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", -108.);
  set_default_double_param("rMin1", 2.2);
  set_default_double_param("rMax1", 65.6);
  set_default_double_param("rMin2", 2.6);
  set_default_double_param("rMax2", 77.5);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 180.);
  set_default_double_param("rot_z", 0.);

  set_default_string_param("material", "G4_PbWO4");
  set_default_string_param("mappingtower", "");
  set_default_string_param("mapping4x4", "");
  return;
}

void PHG4HybridHomogeneousCalorimeterSubsystem::SetTowerMappingFile(const std::string& filename)
{
  set_string_param("mappingtower", filename);
}

void PHG4HybridHomogeneousCalorimeterSubsystem::SetProjectiveGeometry(const std::string& filename1, const std::string& filename2)
{
  set_string_param("mappingtower", filename1);
  set_string_param("mapping4x4", filename2);
  set_int_param("projective", 1);
}
