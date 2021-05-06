#include "SimpleNtuple.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>

#include <sstream>
#include <utility>  // for pair

using namespace std;

SimpleNtuple::SimpleNtuple(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_HistoManager(nullptr)
  , m_Ntup(nullptr)
  , m_Outfile(nullptr)
  , m_Filename(filename)
{
}

SimpleNtuple::~SimpleNtuple()
{
  delete m_HistoManager;
}

int SimpleNtuple::Init(PHCompositeNode *)
{
  m_HistoManager = new Fun4AllHistoManager(Name());
  m_Outfile = new TFile(m_Filename.c_str(), "RECREATE");
  m_Ntup = new TNtuple("hitntup", "G4Hits", "absorber:detid:x0:y0:z0:x1:y1:z1:edep");
  //  ntup->SetDirectory(0);
  TH1 *h1 = new TH1F("edep1GeV", "edep 0-1GeV", 1000, 0, 1);
  m_ElossVec.push_back(h1);
  h1 = new TH1F("edep100GeV", "edep 0-100GeV", 1000, 0, 100);
  m_ElossVec.push_back(h1);
  return 0;
}

int SimpleNtuple::process_event(PHCompositeNode *topNode)
{
  ostringstream nodename;
  set<string>::const_iterator iter;
  vector<TH1 *>::const_iterator eiter;
  for (iter = m_NodePostfixSet.begin(); iter != m_NodePostfixSet.end(); ++iter)
  {
    int detid = (m_DetIdMap.find(*iter))->second;
    nodename.str("");
    nodename << "G4HIT_" << *iter;
    PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (hits)
    {
      double esum = 0;
      PHG4HitContainer::ConstRange hit_range = hits->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)

      {
        esum += hit_iter->second->get_edep();
        m_Ntup->Fill(detid,
                     (int) hit_iter->second->get_layer(), // get_layer is unsigned so cast it back to an int for negative values
                     hit_iter->second->get_x(0),
                     hit_iter->second->get_y(0),
                     hit_iter->second->get_z(0),
                     hit_iter->second->get_x(1),
                     hit_iter->second->get_y(1),
                     hit_iter->second->get_z(1),
                     hit_iter->second->get_edep());
      }
      for (eiter = m_ElossVec.begin(); eiter != m_ElossVec.end(); ++eiter)
      {
        (*eiter)->Fill(esum);
      }
    }
  }
  return 0;
}

int SimpleNtuple::End(PHCompositeNode *topNode)
{
  m_Outfile->cd();
  m_Ntup->Write();
  m_Outfile->Write();
  m_Outfile->Close();
  delete m_Outfile;
  m_HistoManager->dumpHistos(m_Filename, "UPDATE");
  return 0;
}

void SimpleNtuple::AddNode(const std::string &name, const int detid)
{
  m_NodePostfixSet.insert(name);
  m_DetIdMap[name] = detid;
  return;
}
