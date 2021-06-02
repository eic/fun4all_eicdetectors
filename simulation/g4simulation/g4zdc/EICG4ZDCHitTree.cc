// - 1/June/2021 TTree production for ZDC G4 hits     Shima Shimizu
//
#include "EICG4ZDCHitTree.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>           // for SubsysReco

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <sstream>
#include <utility>    

using namespace std;

EICG4ZDCHitTree::EICG4ZDCHitTree(const std::string &name, const std::string &filename) 
  : SubsysReco(name)
  , nblocks(0)
  , hm(nullptr)
  , _filename(filename)
  , tree(nullptr)
  , outfile(nullptr){
}

EICG4ZDCHitTree::~EICG4ZDCHitTree(){
  //  delete ntup;
  delete hm;
}

int EICG4ZDCHitTree::Init(PHCompositeNode *)
 {
   hm = new Fun4AllHistoManager(Name());
   outfile = new TFile(_filename.c_str(), "RECREATE");
   tree = new TTree("zdchit", "Collection of EICG4ZDC G4Hits");

   tree->Branch("Nhit", &Nhit, "Nhit/I");
   tree->Branch("layerType", &layerType);
   tree->Branch("layerID", &layerID);
   tree->Branch("xID", &xID);
   tree->Branch("yID", &yID);
   tree->Branch("x0", &x0);
   tree->Branch("y0", &y0);
   tree->Branch("z0", &z0);
   tree->Branch("x1", &x1);
   tree->Branch("y1", &y1);
   tree->Branch("z1", &z1);
   tree->Branch("time0", &time0);
   tree->Branch("time1", &time1);
   tree->Branch("edep", &edep);

   Nhit=0; 
   layerID=0;
   layerType=0;
   xID=0;
   yID=0;
   x0=0;
   y0=0;
   z0=0;
   x1=0;
   y1=0;
   z1=0;
   time0=0;
   time1=0;
   edep=0;

   return 0;
 }
 
 int EICG4ZDCHitTree::process_event(PHCompositeNode *topNode)
 {
   ostringstream nodename;
   set<string>::const_iterator iter;
   vector<TH1 *>::const_iterator eiter;
   for (iter = _node_postfix.begin(); iter != _node_postfix.end(); ++iter)
   {

     nodename.str("");
     nodename << "G4HIT_" << *iter;
     PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
     if (!hits) return 0;
     
     std::vector<int> v_layerType;
     std::vector<int> v_layerID;
     std::vector<int> v_xID;
     std::vector<int> v_yID;
     std::vector<float> v_x0;
     std::vector<float> v_y0;
     std::vector<float> v_z0;
     std::vector<float> v_x1;
     std::vector<float> v_y1;
     std::vector<float> v_z1;
     std::vector<float> v_time0;
     std::vector<float> v_time1;
     std::vector<float> v_edep;

     Nhit = hits->size();

     PHG4HitContainer::ConstRange hit_range = hits->getHits();
     for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++){
       if(hit_iter->second->get_hit_type()<0) continue;

       v_layerType.push_back(hit_iter->second->get_hit_type());
       v_layerID.push_back(hit_iter->second->get_layer());
       v_xID.push_back(hit_iter->second->get_index_i());
       v_yID.push_back(hit_iter->second->get_index_j());
       v_x0.push_back(hit_iter->second->get_x(0));
       v_y0.push_back(hit_iter->second->get_y(0));
       v_z0.push_back(hit_iter->second->get_z(0));
       v_x1.push_back(hit_iter->second->get_x(1));
       v_y1.push_back(hit_iter->second->get_y(1));
       v_z1.push_back(hit_iter->second->get_z(1));
       v_time0.push_back(hit_iter->second->get_t(0));
       v_time1.push_back(hit_iter->second->get_t(1));
       v_edep.push_back(hit_iter->second->get_edep());
     }

     layerType = &v_layerType;
     layerID = &v_layerID;
     xID = &v_xID;
     yID = &v_yID;
     x0 = &v_x0;
     y0 = &v_y0;
     z0 = &v_z0;
     x1 = &v_x1;
     y1 = &v_y1;
     z1 = &v_z1;
     time0 = &v_time0;
     time1 = &v_time1;
     edep = &v_edep;
     tree->Fill();

     v_layerType.clear();
     v_layerID.clear();
     v_xID.clear();
     v_yID.clear();
     v_x0.clear();
     v_y0.clear();
     v_z0.clear();
     v_x1.clear();
     v_y1.clear();
     v_z1.clear();
     v_time0.clear();
     v_time1.clear();
     v_edep.clear();
   }
   
   return 0;
 }
 
 int EICG4ZDCHitTree::End(PHCompositeNode *topNode)
 {
   outfile->cd();
   tree->Write();
   outfile->Write();
   outfile->Close();
   delete outfile;
   hm->dumpHistos(_filename, "UPDATE");
   return 0;
 }
 
 void EICG4ZDCHitTree::AddNode(const std::string &name, const int detid)
 {
   _node_postfix.insert(name);
   _detid[name] = detid;
   return;
 }
