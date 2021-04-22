#include "AllSiliconTrackerDetector.h"

#include "AllSiliconTrackerDisplayAction.h"
#include "AllSiliconTrackerSubsystem.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Subsystem.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4LogicalVolumeStore.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SolidStore.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <cmath>
#include <iostream>
#include <memory>

class G4VSolid;
class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
AllSiliconTrackerDetector::AllSiliconTrackerDetector(PHG4Subsystem *subsys,
                                                     PHCompositeNode *Node,
                                                     PHParameters *parameters,
                                                     const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<AllSiliconTrackerDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_GDMPath(parameters->get_string_param("GDMPath"))
  , m_Active(m_Params->get_int_param("active"))
  , m_AbsorberActive(parameters->get_int_param("absorberactive"))
{
}

//_______________________________________________________________
int AllSiliconTrackerDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  if (m_Active)
  {
    auto iter = m_ActivePhysicalVolumesSet.find(volume);
    if (iter != m_ActivePhysicalVolumesSet.end())
    {
      return iter->second;
    }
  }
  if (m_AbsorberActive)
  {
    auto iter = m_PassivePhysicalVolumesSet.find(volume);
    if (iter != m_PassivePhysicalVolumesSet.end())
    {
      return iter->second;
    }
  }
  return 0;
}

//_______________________________________________________________
void AllSiliconTrackerDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.SetOverlapCheck(OverlapCheck());
  gdmlParser.Read(m_GDMPath, false);

  // alright the reader just puts everything into G4 Stores - endless fun
  // print the show out, first solids:
  // for (auto i=G4SolidStore::GetInstance()->begin(); i!=G4SolidStore::GetInstance()->end(); ++i)
  //  cout << "solid vol name: " << (*i)->GetName() << endl;
  // for (auto i=G4LogicalVolumeStore::GetInstance()->begin(); i!=G4LogicalVolumeStore::GetInstance()->end(); i++)
  //   cout << "logvol name " << (*i)->GetName() << endl;

  AllSiliconTrackerSubsystem *mysubsys = dynamic_cast<AllSiliconTrackerSubsystem *>(GetMySubsystem());
  for (set<string>::const_iterator its = mysubsys->assembly_iters().first; its != mysubsys->assembly_iters().second; ++its)
  {
    G4AssemblyVolume *avol = reader->GetAssembly(*its);
    if (!avol)
    {
      cout << *its << "not found" << endl;
      continue;
    }
    G4RotationMatrix *rotm = new G4RotationMatrix();
    rotm->rotateX(m_Params->get_double_param("rot_x"));
    rotm->rotateY(m_Params->get_double_param("rot_y"));
    rotm->rotateZ(m_Params->get_double_param("rot_z"));
    G4ThreeVector g4vec(m_Params->get_double_param("place_x"),
                        m_Params->get_double_param("place_y"),
                        m_Params->get_double_param("place_z"));
    avol->MakeImprint(logicWorld, g4vec, rotm, 0, OverlapCheck());
    delete rotm;
    vector<G4VPhysicalVolume *>::iterator it = avol->GetVolumesIterator();
    for (unsigned int i = 0; i < avol->TotalImprintedVolumes(); i++)
    {
      InsertVolumes(*it, insertassemblies);
      ++it;
    }
  }
  for (set<string>::const_iterator its = mysubsys->logvol_iters().first; its != mysubsys->logvol_iters().second; ++its)
  {
    G4LogicalVolume *vol = reader->GetVolume(*its);
    if (!vol)
    {
      cout << *its << "not found" << endl;
      continue;
    }

    G4RotationMatrix *rotm = new G4RotationMatrix();
    rotm->rotateX(m_Params->get_double_param("rot_x"));
    rotm->rotateX(m_Params->get_double_param("rot_y"));
    rotm->rotateX(m_Params->get_double_param("rot_z"));
    G4ThreeVector g4vec(m_Params->get_double_param("place_x"),
                        m_Params->get_double_param("place_y"),
                        m_Params->get_double_param("place_z"));
    G4VPhysicalVolume *phys = new G4PVPlacement(rotm, g4vec,
                                                vol,
                                                G4String(GetName().c_str()),
                                                logicWorld, false, 0, OverlapCheck());
    InsertVolumes(phys, insertlogicalvolumes);
  }
  AddHitNodes(topNode());
  return;
}

void AllSiliconTrackerDetector::InsertVolumes(G4VPhysicalVolume *physvol, const int flag)
{
  static int detid = -9999;
  if (flag == insertassemblies)
  {
    // G4AssemblyVolumes naming convention:
    //     av_WWW_impr_XXX_YYY_ZZZ
    // where:

    //     WWW - assembly volume instance number
    //     XXX - assembly volume imprint number
    //     YYY - the name of the placed logical volume
    //     ZZZ - the logical volume index inside the assembly volume
    // here we enter a new assembly volume and have to set the detector id which
    // stays valid until we hit the next assembly volume
    // the detector id is coded into YYY
    if (physvol->GetName().find("av_") != string::npos && physvol->GetName().find("_impr_") != string::npos)
    {
      detid = -9999;  // reset detid so we see if this is not handled here
      std::vector<std::string> splitname;
      boost::algorithm::split(splitname, physvol->GetName(), boost::is_any_of("_"));
      if (splitname[4].find("AluStrips") != string::npos)
      {
        detid = 100;
      }
      else
      {
        string detprefix[] = {"VstStave", "FstContainerVolume", "BstContainerVolume", "Beampipe"};
        //start layer count at 10 for VstStave
        // at 20 for Fst
        // at 30 for Bst
        // at 40 for beampipe
        int increase = 10;
        for (auto toerase : detprefix)
        {
          size_t pos = splitname[4].find(toerase);
          if (pos != string::npos)
          {
            detid = boost::lexical_cast<int>(splitname[4].erase(pos, toerase.length())) + increase;
            break;
          }
          increase += 10;
        }
      }
    }
    if (detid < 0)
    {
      cout << "AllSiliconTrackerDetector::InsertVolumes detid is " << detid << " vol name " << physvol->GetName() << endl;
      gSystem->Exit(1);
    }
  }
  else
  {
    detid = 1;
  }
  G4LogicalVolume *logvol = physvol->GetLogicalVolume();
  m_DisplayAction->AddLogicalVolume(logvol);
  //    cout << "mat: " << logvol->GetMaterial()->GetName() << endl;
  //  cout << "Adding " << physvol->GetName() << endl;
  if (physvol->GetName().find("MimosaCore") != string::npos)
  {
    m_ActivePhysicalVolumesSet.insert(make_pair(physvol, detid));
  }
  else
  {
    if (m_AbsorberActive)
    {
      m_PassivePhysicalVolumesSet.insert(make_pair(physvol, -detid));
    }
  }
  // G4 10.06 returns unsigned int for GetNoDaughters()
  // lower version return int, need to cast to avoid compiler error
  for (int i = 0; i < (int) logvol->GetNoDaughters(); ++i)
  {
    G4VPhysicalVolume *physvol = logvol->GetDaughter(i);
    // here we decide which volumes are active
    InsertVolumes(physvol, flag);
  }
  return;
}

//_______________________________________________________________
void AllSiliconTrackerDetector::Print(const std::string &what) const
{
  cout << "AllSiliconTracker Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Version 0.1" << endl;
    cout << "Parameters:" << endl;
    m_Params->Print();
  }
  return;
}

int AllSiliconTrackerDetector::get_detid(const G4VPhysicalVolume *physvol, const int whichactive)
{
  if (whichactive > 0)
  {
    return m_ActivePhysicalVolumesSet.find(physvol)->second;
  }
  else if (whichactive < 0)
  {
    return m_PassivePhysicalVolumesSet.find(physvol)->second;
  }
  else
  {
    cout << "AllSiliconTrackerDetector::get_detid invalid whichactive flag = 0" << endl;
    gSystem->Exit(1);
  }
  return 0;
}

void AllSiliconTrackerDetector::AddHitNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  PHNodeIterator dstIter(dstNode);
  if (m_Params->get_int_param("active"))
  {
    map<string, int> nodes;
    string myname;
    if (SuperDetector() != "NONE")
    {
      myname = SuperDetector();
    }
    else
    {
      myname = GetName();
    }
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", myname));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(myname);
      dstNode->addNode(DetNode);
    }
    set<int> DetidSet;
    for (auto iter = m_ActivePhysicalVolumesSet.begin(); iter != m_ActivePhysicalVolumesSet.end(); ++iter)
    {
      DetidSet.insert(iter->second);
    }
    ostringstream g4hitnodeactive;
    for (auto iter = DetidSet.begin(); iter != DetidSet.end(); ++iter)
    {
      g4hitnodeactive.str("");
      string region;
      if (*iter < 20)
      {
        region = "_CENTRAL_";
      }
      else if (*iter < 30)
      {
        region = "_FORWARD_";
      }
      else
      {
        region = "_BACKWARD_";
      }
      g4hitnodeactive << "G4HIT_" << myname << region << *iter;
      nodes.insert(make_pair(g4hitnodeactive.str(), *iter));
    }
    if (m_Params->get_int_param("absorberactive"))
    {
      string g4hitnodename = "G4HIT_ABSORBER_" + myname;
      nodes.insert(make_pair(g4hitnodename, -1));
    }
    for (auto nodename : nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(DetNode, nodename.first);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename.first);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename.first, "PHObject"));
      }
      m_HitContainerMap.insert(make_pair(nodename.second, g4_hits));
    }
  }
  return;
}

PHG4HitContainer *AllSiliconTrackerDetector::get_hitcontainer(const int i)
{
  map<int, PHG4HitContainer *>::const_iterator iter = m_HitContainerMap.find(i);
  if (iter != m_HitContainerMap.end())
  {
    return iter->second;
  }
  cout << "bad index: " << i << endl;
  gSystem->Exit(1);
  exit(1);
}
