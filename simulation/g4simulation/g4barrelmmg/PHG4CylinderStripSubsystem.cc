#include "PHG4CylinderStripSubsystem.h"
#include "PHG4CylinderStripDetector.h"
#include "PHG4CylinderStripSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>    // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>   // for PHG4SteppingAction
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>          // for PHIODataNode
#include <phool/PHNode.h>                // for PHNode
#include <phool/PHNodeIterator.h>        // for PHNodeIterator
#include <phool/PHObject.h>              // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4Types.hh>             // for G4double

#include <cmath>                        // for NAN
#include <iostream>                      // for operator<<, basic_ostream, endl
#include <sstream>

class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4CylinderStripSubsystem::PHG4CylinderStripSubsystem(const std::string &na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
{
  InitializeParameters();
}

PHG4CylinderStripSubsystem::~PHG4CylinderStripSubsystem()
{
}

//_______________________________________________________________________
int PHG4CylinderStripSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  // create detector
  m_Detector = new PHG4CylinderStripDetector(this, topNode, GetParams(), Name(), GetLayer());
  //G4double detlength = GetParams()->get_double_param("length");
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));

    ostringstream nodename;
    ostringstream geonode;
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
      PHNodeIterator iter_run(runNode);
      superSubNode = dynamic_cast<PHCompositeNode *>(iter_run.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        runNode->addNode(superSubNode);
      }
      runNode = superSubNode;

      nodename << "G4HIT_" << SuperDetector();
      geonode << "CYLINDERStripGEOM_" << SuperDetector();
    }

    else
    {
      nodename << "G4HIT_" << Name();
      geonode << "CYLINDERStripGEOM_" << Name();
    }
    PHG4HitContainer *cylinder_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!cylinder_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(cylinder_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }
    cylinder_hits->AddLayer(GetLayer());
    //PHG4CylinderStripGeomContainer *geo = findNode::getClass<PHG4CylinderStripGeomContainer>(topNode, geonode.str());
    //if (!geo)
    //{
    //  geo = new PHG4CylinderStripGeomContainer();
    //  PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode.str(), "PHObject");
    //  runNode->addNode(newNode);
    //}
    //PHG4CylinderStripGeom *mygeom = new PHG4CylinderStripGeomv1(GetParams()->get_double_param("radius"), GetParams()->get_double_param("place_z") - detlength / 2., GetParams()->get_double_param("place_z") + detlength / 2., GetParams()->get_double_param("thickness"));
    //geo->AddLayerGeom(GetLayer(), mygeom);
    m_SteppingAction = new PHG4CylinderStripSteppingAction(m_Detector, GetParams());
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new PHG4CylinderStripSteppingAction(m_Detector, GetParams());
  }
  return 0;
}

//_______________________________________________________________________
int PHG4CylinderStripSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4CylinderStripSubsystem::SetDefaultParameters()
{
  set_default_double_param("length", 100);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("radius", 100);
  set_default_double_param("steplimits", NAN);
  set_default_double_param("gap", 20);
  set_default_double_param("gas1thickness", 0.0020000);
  set_default_double_param("gas2thickness", 0.3000000);
  set_default_double_param("phi0", 0.);
  set_default_double_param("tmin", NAN);
  set_default_double_param("tmax", NAN);
  set_default_double_param("deadzone", 0);

  set_default_int_param("lengthviarapidity", 0);
  set_default_int_param("lightyield", 0);
  set_default_int_param("use_g4steps", 0);
  set_default_int_param("use_2Dreadout", 0);
  set_default_int_param("nhit", 1);

  set_default_string_param("gas", "G4_Ar");
}

PHG4Detector *
PHG4CylinderStripSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4CylinderStripSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  if (!BeginRunExecuted())
  {
    cout << "Need to execute BeginRun() before parameter printout is meaningful" << endl;
    cout << "To do so either run one or more events or on the command line execute: " << endl;
    cout << "Fun4AllServer *se = Fun4AllServer::instance();" << endl;
    cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << endl;
    cout << "g4->InitRun(se->topNode());" << endl;
    cout << "PHG4CylinderStripSubsystem *cyl = (PHG4CylinderStripSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << endl;
    cout << "cyl->Print()" << endl;
    return;
  }
  GetParams()->Print();
  if (m_SteppingAction)
  {
    m_SteppingAction->Print(what);
  }
  return;
}
