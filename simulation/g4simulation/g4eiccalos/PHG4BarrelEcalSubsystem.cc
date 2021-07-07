#include "PHG4BarrelEcalSubsystem.h"
#include "PHG4BarrelEcalDetector.h"
#include "PHG4BarrelEcalDetector.h"
#include "PHG4BarrelEcalDisplayAction.h"
#include "PHG4BarrelEcalSteppingAction.h"

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

#include <TSystem.h>

#include <cstdlib>  // for getenv
#include <set>      // for set
#include <sstream>

class PHG4Detector;
 
//_______________________________________________________________________
PHG4BarrelEcalSubsystem::PHG4BarrelEcalSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4BarrelEcalSubsystem::~PHG4BarrelEcalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4BarrelEcalSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4BarrelEcalDisplayAction(Name());
  // create detector
  m_Detector = new PHG4BarrelEcalDetector(this, topNode, GetParams(), Name());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());
  
  std::set<std::string> nodes;

  if (GetParams()->get_int_param("active"))
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
    // create hit output nodes
    std::string detector_suffix = SuperDetector();
    if (detector_suffix == "NONE")
    {
      detector_suffix = Name();
    }
    m_HitNodeName = "G4HIT_" + detector_suffix;
    nodes.insert(m_HitNodeName);
    m_AbsorberNodeName = "G4HIT_ABSORBER_" + detector_suffix;
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodes.insert(m_AbsorberNodeName);
    }
    m_SupportNodeName = "G4HIT_SUPPORT_" + detector_suffix;
    if (GetParams()->get_int_param("supportactive"))
    {
      nodes.insert(m_SupportNodeName);
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
    PHG4BarrelEcalSteppingAction* tmp = new PHG4BarrelEcalSteppingAction(m_Detector, GetParams());
    tmp->SetHitNodeName(m_HitNodeName);
    tmp->SetAbsorberNodeName(m_AbsorberNodeName);
    tmp->SetSupportNodeName(m_SupportNodeName);
    m_SteppingAction = tmp;
  }
  return 0;
}

//_______________________________________________________________________
int PHG4BarrelEcalSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4BarrelEcalSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4BarrelEcalSubsystem::SetDefaultParameters()
{
  std::ostringstream mappingfilename;
  const char* calibroot = getenv("CALIBRATIONROOT");
  if (calibroot)
  {
    mappingfilename << calibroot;
  }
  else
  {
    std::cout << "no CALIBRATIONROOT environment variable" << std::endl;
    gSystem->Exit(1);
  }
  mappingfilename << "BarrelEcal/mapping/towerMap_BEMC_v001.txt";
  set_default_string_param("mapping_file", mappingfilename.str());
  set_default_double_param("radius", 85.);
  set_default_double_param("Length", 298.94);
  set_default_double_param("max_radius", 138.);


  return;
}

void PHG4BarrelEcalSubsystem::SetTowerMappingFile(const std::string& filename)
{
  set_string_param("mapping_file", filename);
  
}

void PHG4BarrelEcalSubsystem::Print(const std::string& what) const
{
  std::cout << Name() << " Parameters: " << std::endl;
  if (!BeginRunExecuted())
  {
    std::cout << "Need to execute BeginRun() before parameter printout is meaningful" << std::endl;
    std::cout << "To do so either run one or more events or on the command line execute: " << std::endl;
    std::cout << "Fun4AllServer *se = Fun4AllServer::instance();" << std::endl;
    std::cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << std::endl;
    std::cout << "g4->InitRun(se->topNode());" << std::endl;
    std::cout << "PHG4BarrelEcalSubsystem *fhcal = (PHG4BarrelEcalSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << std::endl;
    std::cout << "fhcal->Print()" << std::endl;
    return;
  }
  GetParams()->Print();
  if (m_SteppingAction)
  {
    m_SteppingAction->Print(what);
  }
  return;
}