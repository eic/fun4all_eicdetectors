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

class EICPIDParticleContainer : public PHObject
{
 public:
  typedef std::map<EICPIDDefs::keytype, EICPIDParticle*> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<unsigned int>::const_iterator LayerIter;

  EICPIDParticleContainer();

  ~EICPIDParticleContainer() override {}

  void Reset() override;

  void identify(std::ostream& os = std::cout) const override;

  ConstIterator AddPIDParticle(EICPIDParticle* newhit);

  Iterator findOrAddPIDParticle(EICPIDDefs::keytype key);

  EICPIDParticle* findEICPIDParticle(EICPIDDefs::keytype key);

  //! return all hist
  ConstRange getPIDParticles(void) const;

  unsigned int size(void) const
  {
    return m_particleMap.size();
  }

 protected:
  Map m_particleMap;

  ClassDefOverride(EICPIDParticleContainer, 1)
};

#endif
