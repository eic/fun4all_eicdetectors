// - 1/June/2021 TTree production for B0 G4 hits     Shima Shimizu
//
#include "EICG4B0HitTree.h"

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

EICG4B0HitTree::EICG4B0HitTree(const std::string &name, const std::string &filename)
    : SubsysReco(name)
    , _filename(filename)
    , tree(nullptr)
    , outfile(nullptr)
    , Nhit(0)
{
}

EICG4B0HitTree::~EICG4B0HitTree(){
}

int EICG4B0HitTree::Init(PHCompositeNode *)
{
    outfile = new TFile(_filename.c_str(), "RECREATE");
    tree = new TTree("b0hit", "Collection of EICG4 B0 G4Hits");

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

    return 0;
}

int EICG4B0HitTree::process_event(PHCompositeNode *topNode)
{
    ostringstream nodename;
    set<string>::const_iterator iter;
    for (iter = _node_postfix.begin(); iter != _node_postfix.end(); ++iter)
    {
        nodename.str("");
        nodename << "G4HIT_" << *iter;
        PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
        if (!hits) return 0;

        Nhit = hits->size();

        PHG4HitContainer::ConstRange hit_range = hits->getHits();
        for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++){
//            cout << " ---> !!! hit_iter->first = " << hit_iter->first << " hit_iter->second->get_hit_type() = " << hit_iter->second->get_hit_type() << endl;
//            hit_iter->second->identify();
//            hit_iter->second->print();
            if(hit_iter->second->get_hit_type()<0) continue;

            layerType.push_back(hit_iter->second->get_hit_type());
            layerID.push_back(hit_iter->second->get_layer());
            xID.push_back(hit_iter->second->get_index_i());
            yID.push_back(hit_iter->second->get_index_j());
            x0.push_back(hit_iter->second->get_x(0));
            y0.push_back(hit_iter->second->get_y(0));
            z0.push_back(hit_iter->second->get_z(0));
            x1.push_back(hit_iter->second->get_x(1));
            y1.push_back(hit_iter->second->get_y(1));
            z1.push_back(hit_iter->second->get_z(1));
            time0.push_back(hit_iter->second->get_t(0));
            time1.push_back(hit_iter->second->get_t(1));
            edep.push_back(hit_iter->second->get_edep());
        }

        tree->Fill();

        layerType.clear();
        layerID.clear();
        xID.clear();
        yID.clear();
        x0.clear();
        y0.clear();
        z0.clear();
        x1.clear();
        y1.clear();
        z1.clear();
        time0.clear();
        time1.clear();
        edep.clear();
    }
    return 0;
}

int EICG4B0HitTree::End(PHCompositeNode *topNode)
{
    outfile->cd();
    tree->Write();
    outfile->Write();
    outfile->Close();
    delete outfile;
    return 0;
}

void EICG4B0HitTree::AddNode(const std::string &name, const int detid)
{
    _node_postfix.insert(name);
    _detid[name] = detid;
    return;
}

