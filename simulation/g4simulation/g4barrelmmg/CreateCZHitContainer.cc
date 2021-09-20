#include "CreateCZHitContainer.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>           // for SubsysReco
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>          // for PHIODataNode
#include <phool/PHNode.h>                // for PHNode
#include <phool/PHNodeIterator.h>        // for PHNodeIterator
#include <phool/PHObject.h>              // for PHObject
#include <phool/getClass.h>
#include <TSystem.h>
#include <TMath.h>

#include <sstream>
#include <utility>                        // for pair

using namespace std;

CreateCZHitContainer::CreateCZHitContainer(const std::string &name)
  : SubsysReco("CreateCZHitContainer"+name)
  , _node_postfix(name)
{
}

CreateCZHitContainer::~CreateCZHitContainer()
{
}

int CreateCZHitContainer::InitRun(PHCompositeNode *topNode)
{
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth_container)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  ostringstream nodename;
  nodename.str("");
  nodename << "G4HIT_CZ" << _node_postfix;
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", _node_postfix));
  PHG4HitContainer *cylinder_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
  if (!cylinder_hits)
  {
    dstNode->addNode(new PHIODataNode<PHObject>(cylinder_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
  }
  //cylinder_hits->AddLayer(GetLayer());
  return Fun4AllReturnCodes::EVENT_OK;
}

int CreateCZHitContainer::process_event(PHCompositeNode *topNode)
{
  ostringstream nodename;
  nodename.str("");
  nodename << "G4HIT_" << _node_postfix;
  PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  nodename.str("");
  nodename.clear();
  nodename << "G4HIT_CZ" << _node_postfix;
  PHG4HitContainer *cylinder_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  if (!cylinder_hits){
	cerr << "can't find " << nodename.str() << endl;
  	gSystem->Exit(1);
  }
  if (hits)
  {
	PHG4TruthInfoContainer::ConstRange itr_range = _truth_container->GetParticleRange();
	for (PHG4TruthInfoContainer::ConstIterator itr = itr_range.first;
       itr != itr_range.second; ++itr)
	{
	  PHG4Particle* particle = itr->second;
	  if (!particle) continue;
	  int trkid = particle->get_track_id();
      for (PHG4HitContainer::LayerIter layerit =
                 hits->getLayers().first;
             layerit != hits->getLayers().second; layerit++){
        cylinder_hits->AddLayer(*layerit);
        PHG4HitContainer::ConstRange hit_range = hits->getHits(*layerit);
        _hit_CZ = nullptr;
        _hit_C = nullptr;
        _hit_Z = nullptr;
        for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
        {
	      if (hit_iter->second->get_trkid() != trkid) continue;
	      if (hit_iter->second->get_hit_type() == 1)
	        _hit_C = hit_iter->second;
		  else if (hit_iter->second->get_hit_type() == 2)
			_hit_Z = hit_iter->second;
        }
		_hit_CZ = merge_hits (_hit_C, _hit_Z);
		if (_hit_CZ){
	      cylinder_hits -> AddHit(*layerit, _hit_CZ);
		}
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

PHG4Hit* CreateCZHitContainer::merge_hits(PHG4Hit* h1, PHG4Hit* h2){
  if (h1==nullptr || h2==nullptr) return nullptr;
  if (h1->get_trkid() != h2->get_trkid()) return nullptr;
  if (h1->get_layer() != h2->get_layer()) return nullptr;
  if (!(h1->get_hit_type() == 1 && h2->get_hit_type() == 2)) return nullptr;
  if (!(h1->get_eion() > 0 && h2->get_eion() > 0)) return nullptr;
  PHG4Hit* merged = new PHG4Hitv1();
  merged->set_layer(h1->get_layer());

  float R_in = TMath::Sqrt(h1->get_x(0)*h1->get_x(0) + h1->get_y(0)*h1->get_y(0));
  float R_out = TMath::Sqrt(h2->get_x(1)*h2->get_x(1) + h2->get_y(1)*h2->get_y(1));
  float phi_measured = TMath::ATan2(h2->get_avg_y(), h2->get_avg_x());

  merged->set_x(0, R_in * TMath::Cos(phi_measured));
  merged->set_y(0, R_in * TMath::Sin(phi_measured));
  merged->set_z(0, h1->get_z(0));
  merged->set_x(1, R_out * TMath::Cos(phi_measured));
  merged->set_y(1, R_out * TMath::Sin(phi_measured));
  merged->set_z(1, h1->get_z(1));
  float t0 = h1->get_t(0);
  if (h2->get_t(0) < t0) t0=h2->get_t(0);
  float t1 = h1->get_t(1);
  if (h2->get_t(1) > t1) t1=h2->get_t(1);
  merged->set_t(0, t0);
  merged->set_t(1, t1);
  merged->set_trkid(h1->get_trkid());
  return merged;
}
