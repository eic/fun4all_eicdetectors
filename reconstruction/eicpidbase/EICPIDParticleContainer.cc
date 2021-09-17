#include "EICPIDParticleContainer.h"

#include <phool/phool.h>

#include <TSystem.h>

#include <cstdlib>
#include "EICPIDParticle.h"
#include "EICPIDParticlev1.h"

using namespace std;

EICPIDParticleContainer::EICPIDParticleContainer()
  : id(-1), hitmap(), layers()
{
}

EICPIDParticleContainer::EICPIDParticleContainer(const std::string &nodename)
  : id(EICPIDDefs::get_volume_id(nodename)), hitmap(), layers()
{
}

void
EICPIDParticleContainer::Reset()
{
   while(hitmap.begin() != hitmap.end())
     {
       delete hitmap.begin()->second;
       hitmap.erase(hitmap.begin());
     }
  return;
}

void
EICPIDParticleContainer::identify(ostream& os) const
{
   ConstIterator iter;
   os << "Number of hits: " << size() << endl;
   for (iter = hitmap.begin(); iter != hitmap.end(); ++iter)
     {
       os << "hit key 0x" << hex << iter->first << dec << endl;
       (iter->second)->identify();
     }
   set<unsigned int>::const_iterator siter;
   os << "Number of layers: " << num_layers() << endl;
   for (siter = layers.begin(); siter != layers.end(); ++siter)
     {
       os << "layer : " << *siter << endl;
     }
  return;
}

EICPIDDefs::keytype
EICPIDParticleContainer::getmaxkey(const unsigned int detid)
{
  ConstRange miter = getHits(detid);
  // first handle 2 special cases where there is no hit in the current layer
  // no hits in this layer and higher layers (lower layers can contain hits)
    if (miter.first == hitmap.end())
    {
      return 0;
    }
  // no hits in this layer - but hits in higher layers
  if (miter.first == miter.second)
    {
      return  0;
    }
  EICPIDDefs::keytype detidlong = detid;
  EICPIDDefs::keytype shiftval = detidlong << EICPIDDefs::hit_idbits;
  ConstIterator lastlayerentry = miter.second;
  --lastlayerentry;
  EICPIDDefs::keytype iret = lastlayerentry->first - shiftval; // subtract layer mask
  return iret;
}


EICPIDDefs::keytype
EICPIDParticleContainer::genkey(const unsigned int detid)
{
  EICPIDDefs::keytype detidlong = detid;
  if ((detidlong >> EICPIDDefs::keybits) > 0)
    {
      cout << PHWHERE << " detector id too large: " << detid << endl;
      gSystem->Exit(1);
    }
  EICPIDDefs::keytype shiftval = detidlong << EICPIDDefs::hit_idbits;
  //  cout << "max index: " << (detminmax->second)->first << endl;
  // after removing hits with no energy deposition, we have holes
  // in our hit ranges. This construct will get us the last hit in
  // a layer and return it's hit id. Adding 1 will put us at the end of this layer
  EICPIDDefs::keytype hitid = getmaxkey(detid);
  hitid++;
  EICPIDDefs::keytype newkey = hitid | shiftval;
  if (hitmap.find(newkey) != hitmap.end())
    {
      cout << PHWHERE << " duplicate key: 0x" 
           << hex << newkey << dec 
	   << " for detector " << detid 
	   << " hitmap.size: " << hitmap.size()
	   << " hitid: " << hitid << " exiting now" << endl;
            exit(1);
    }
  return newkey;
}

EICPIDParticleContainer::ConstIterator
EICPIDParticleContainer::AddHit(EICPIDParticle *newhit)
{
  EICPIDDefs::keytype key = newhit->get_hit_id();
  if (hitmap.find(key) != hitmap.end())
    {
      cout << "hit with id  0x" << hex << key << dec << " exists already" << endl;
      return hitmap.find(key);
    }
  EICPIDDefs::keytype detidlong = key >>  EICPIDDefs::hit_idbits;
  unsigned int detid = detidlong;
  layers.insert(detid);
  return hitmap.insert( std::make_pair( key, newhit ) ).first;
}

EICPIDParticleContainer::ConstIterator
EICPIDParticleContainer::AddHit(const unsigned int detid, EICPIDParticle *newhit)
{
  EICPIDDefs::keytype key = genkey(detid);
  layers.insert(detid);
  newhit->set_hit_id(key);
  return hitmap.insert( std::make_pair( key, newhit ) ).first;
}

EICPIDParticleContainer::ConstRange EICPIDParticleContainer::getHits(const unsigned int detid) const
{
  EICPIDDefs::keytype detidlong = detid;
  if ((detidlong >> EICPIDDefs::keybits) > 0)
    {
      cout << " detector id too large: " << detid << endl;
      exit(1);
    }
  EICPIDDefs::keytype keylow = detidlong << EICPIDDefs::hit_idbits;
  EICPIDDefs::keytype keyup = ((detidlong + 1) << EICPIDDefs::hit_idbits) -1 ;
  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylow);
  retpair.second = hitmap.upper_bound(keyup);
  return retpair;
}

EICPIDParticleContainer::ConstRange EICPIDParticleContainer::getHits( void ) const
{ return std::make_pair( hitmap.begin(), hitmap.end() ); }


EICPIDParticleContainer::Iterator EICPIDParticleContainer::findOrAddHit(EICPIDDefs::keytype key)
{
  EICPIDParticleContainer::Iterator it = hitmap.find(key);
  if(it == hitmap.end())
  {
    hitmap[key] = new EICPIDParticlev1();
    it = hitmap.find(key);
    EICPIDParticle* mhit = it->second;
    mhit->set_hit_id(key);
    mhit->set_edep(0.);
    layers.insert(mhit->get_layer()); // add layer to our set of layers
  }
  return it;
}

EICPIDParticle* EICPIDParticleContainer::findHit(EICPIDDefs::keytype key)
{
  EICPIDParticleContainer::ConstIterator it = hitmap.find(key);
  if(it != hitmap.end())
    {
      return it->second;
    }
    
  return nullptr;
}

void
EICPIDParticleContainer::RemoveZeroEDep()
{
  //  unsigned int hitsbef = hitmap.size();
  Iterator itr = hitmap.begin();
  Iterator last = hitmap.end();
  for (; itr != last; )
    {
      EICPIDParticle *hit = itr->second;
      if (hit->get_edep() == 0)
        {
          delete hit;
          hitmap.erase(itr++);
        }
      else
        {
          ++itr;
        }
    }
//   unsigned int hitsafter = hitmap.size();
//   cout << "hist before: " << hitsbef
//        << ", hits after: " << hitsafter << endl;
  return;
}

