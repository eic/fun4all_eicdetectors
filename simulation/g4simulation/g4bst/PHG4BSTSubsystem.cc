#include "PHG4BSTSubsystem.h"
#include "PHG4BSTDetector.h"
#include "PHG4BSTDisplayAction.h"
#include "PHG4BSTSteppingAction.h"

#include <g4main/PHG4DisplayAction.h>       // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>      // for PHG4SteppingAction
#include <g4main/PHG4Subsystem.h>           // for PHG4Subsystem
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>             // for PHIODataNode
#include <phool/PHNode.h>                   // for PHNode
#include <phool/PHNodeIterator.h>           // for PHNodeIterator
#include <phool/PHObject.h>                 // for PHObject
#include <phool/getClass.h>

#include <set>                              // for set
#include <sstream>

class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4BSTSubsystem::PHG4BSTSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
  , active(1)
  , absorber_active(0)
  , blackhole(0)
  , detector_type(name)
  , mappingfile_("")
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4BSTSubsystem::~PHG4BSTSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4BSTSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4BSTDisplayAction(Name());
  // create detector
  m_Detector = new PHG4BSTDetector(this, topNode, GetParams(), Name());
  m_Detector->SetActive(active);
  m_Detector->SetAbsorberActive(absorber_active);
  m_Detector->BlackHole(blackhole);
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());
  // m_Detector->SetTowerMappingFile(mappingfile_);
  m_Detector->SuperDetector(SuperDetector());





  std::set<std::string> nodes;
  if (active)
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* DetNode = dstNode;
    if (SuperDetector() != "NONE")
    {
      DetNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(DetNode);
      }
    }
    for(int ilay=0; ilay<6; ilay++)
    {
      // create hit output node
      std::string nodename;
      if (SuperDetector() != "NONE")
      {
        nodename = "G4HIT_" + SuperDetector() + "_" + std::to_string(ilay);
      }
      else
      {
        nodename = "G4HIT_" + Name() + "_" + std::to_string(ilay);
      }
      nodes.insert(nodename);
    }

    if (absorber_active)
    {
      std::string nodename;
      if (SuperDetector() != "NONE")
      {
        nodename = "G4HIT_ABSORBER_" + SuperDetector();
      }
      else
      {
        nodename = "G4HIT_ABSORBER_" + Name();
      }
      nodes.insert(nodename);
    }

    for (auto thisnode : nodes)
    {
      PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(topNode, thisnode);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(thisnode);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, thisnode, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4BSTSteppingAction(m_Detector,absorber_active);
    m_Detector->SetSteppingAction(dynamic_cast<PHG4BSTSteppingAction*>(m_SteppingAction));
  }

  return 0;
}

//_______________________________________________________________________
int PHG4BSTSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4BSTSubsystem::GetDetector() const
{
  return m_Detector;
}


void PHG4BSTSubsystem::SetDefaultParameters()
{
  set_default_double_param("layer_backing_thickness", 0.0);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 375.);
  set_default_double_param("tower_dx", 0.3);
  set_default_double_param("tower_dy", 0.3);
  set_default_double_param("tower_dz", 150.);
  set_default_double_param("dz", 150.);
  set_default_double_param("rMin1", 20.);
  set_default_double_param("rMax1", 220.);
  set_default_double_param("rMin2", 20.);
  set_default_double_param("rMax2", 220.);
  set_default_int_param("do_internal_supports", 1);
  set_default_int_param("do_external_supports", 1);
  set_default_int_param("use_bent_wafer_sagittas_default", 0);
  set_default_int_param("use_bent_wafer_sagittas_mod", 0);
  set_default_int_param("use_ECCE_with_OuterStave", 0);
  set_default_int_param("use_EPIC_setup", 0);
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


void PHG4BSTSubsystem::SetTowerMappingFile(const std::string& filename)
{
  // set_string_param("mapping_file", filename);
  // set_string_param("mapping_file_md5", PHG4Utils::md5sum(get_string_param("mapping_file")));
  mappingfile_ = filename;
}