// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICPID_EICPIDParticleCONTAINER_H
#define EICPID_EICPIDParticleCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include "EICPIDDefs.h"

class EICPIDParticle;

class EICPIDParticleContainer: public PHObject
{

  public:
  typedef std::map<EICPIDDefs::keytype, EICPIDParticle *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<unsigned int>::const_iterator LayerIter;

  EICPIDParticleContainer(); //< used only by ROOT for DST readback
  EICPIDParticleContainer(const std::string &nodename);

  ~EICPIDParticleContainer() override {}

  void Reset() override;

  void identify(std::ostream& os = std::cout) const override;

  //! container ID should follow definition of EICPIDDefs::get_volume_id(DST nodename)
  void SetID(int i) {id = i;}
  int GetID() const {return id;}
  
  ConstIterator AddHit(EICPIDParticle *newhit);

  ConstIterator AddHit(const unsigned int detid, EICPIDParticle *newhit);
  
  Iterator findOrAddHit(EICPIDDefs::keytype key);

  EICPIDParticle* findHit(EICPIDDefs::keytype key );

  EICPIDDefs::keytype genkey(const unsigned int detid);

  //! return all hits matching a given detid
  ConstRange getHits(const unsigned int detid) const;

  //! return all hist
  ConstRange getHits( void ) const;

  unsigned int size( void ) const
  { return hitmap.size(); }
  unsigned int num_layers(void) const
  { return layers.size(); }
  std::pair<LayerIter, LayerIter> getLayers() const
     { return make_pair(layers.begin(), layers.end());} 
  void AddLayer(const unsigned int ilayer) {layers.insert(ilayer);}
  void RemoveZeroEDep();
  EICPIDDefs::keytype getmaxkey(const unsigned int detid);

 protected:

  int id; //< unique identifier from hash of node name. Defined following EICPIDDefs::get_volume_id
  Map hitmap;
  std::set<unsigned int> layers; // layers is not reset since layers must not change event by event

  ClassDefOverride(EICPIDParticleContainer,1)
};

#endif
