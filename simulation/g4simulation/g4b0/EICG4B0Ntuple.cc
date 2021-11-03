//  Ntuple production for B0 G4 hits     Sasha Bylinkin
//
#include "EICG4B0Ntuple.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>

#include <sstream>
#include <utility>

EICG4B0Ntuple::EICG4B0Ntuple(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , hm(nullptr)
  , _filename(filename)
  , ntup(nullptr)
  , outfile(nullptr)
{
}

EICG4B0Ntuple::~EICG4B0Ntuple()
{
  //  delete ntup;
  delete hm;
}

int EICG4B0Ntuple::Init(PHCompositeNode *)
{
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(_filename.c_str(), "RECREATE");
  ntup = new TNtuple("b0hit", "B0s", "layer:x0:y0:z0:x1:y1:z1:t0:t1:edep");
  //  ntup->SetDirectory(0);
  return 0;
}

int EICG4B0Ntuple::process_event(PHCompositeNode *topNode)
{
  std::ostringstream nodename;
  std::set<std::string>::const_iterator iter;
  for (iter = _node_postfix.begin(); iter != _node_postfix.end(); ++iter)
  {
    //int detid = (_detid.find(*iter))->second;
    nodename.str("");
    nodename << "G4HIT_" << *iter;
    PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!hits) return 0;

    double esum = 0;
    //          double numhits = hits->size();
    //          nhits[i]->Fill(numhits);
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
      if (hit_iter->second->get_hit_type() < 0) continue;

      esum += hit_iter->second->get_edep();
      ntup->Fill(hit_iter->second->get_layer(),
                 hit_iter->second->get_x(0),
                 hit_iter->second->get_y(0),
                 hit_iter->second->get_z(0),
                 hit_iter->second->get_x(1),
                 hit_iter->second->get_y(1),
                 hit_iter->second->get_z(1),
                 hit_iter->second->get_t(0),
                 hit_iter->second->get_t(1),
                 hit_iter->second->get_edep());
    }
  }

  return 0;
}

int EICG4B0Ntuple::End(PHCompositeNode *topNode)
{
  outfile->cd();
  ntup->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(_filename, "UPDATE");
  return 0;
}

void EICG4B0Ntuple::AddNode(const std::string &name, const int detid)
{
  _node_postfix.insert(name);
  _detid[name] = detid;
  return;
}
