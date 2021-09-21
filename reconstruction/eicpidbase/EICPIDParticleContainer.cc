#include "EICPIDParticleContainer.h"

#include <phool/phool.h>

#include <TSystem.h>

#include <cstdlib>
#include "EICPIDParticle.h"
#include "EICPIDParticlev1.h"

using namespace std;

EICPIDParticleContainer::EICPIDParticleContainer()
{
}

void EICPIDParticleContainer::Reset()
{
  while (m_particleMap.begin() != m_particleMap.end())
  {
    delete m_particleMap.begin()->second;
    m_particleMap.erase(m_particleMap.begin());
  }
  return;
}

void EICPIDParticleContainer::identify(ostream& os) const
{
  ConstIterator iter;
  os << "Number of PIDParticles: " << size() << endl;
  for (iter = m_particleMap.begin(); iter != m_particleMap.end(); ++iter)
  {
    os << "PIDParticles ID " << iter->first << ": ";
    (iter->second)->identify();
  }
  return;
}

EICPIDParticleContainer::ConstIterator
EICPIDParticleContainer::AddPIDParticle(EICPIDParticle* newhit)
{
  EICPIDDefs::keytype key = newhit->get_id();
  if (m_particleMap.find(key) != m_particleMap.end())
  {
    cout << "hit with id " << key << " exists already" << endl;
    return m_particleMap.find(key);
  }
  return m_particleMap.insert(std::make_pair(key, newhit)).first;
}

EICPIDParticleContainer::ConstRange EICPIDParticleContainer::getPIDParticles(void) const
{
  return std::make_pair(m_particleMap.begin(), m_particleMap.end());
}

EICPIDParticleContainer::Iterator EICPIDParticleContainer::findOrAddPIDParticle(EICPIDDefs::keytype key)
{
  EICPIDParticleContainer::Iterator it = m_particleMap.find(key);
  if (it == m_particleMap.end())
  {
    m_particleMap[key] = new EICPIDParticlev1();
    it = m_particleMap.find(key);
    EICPIDParticle* mhit = it->second;
    mhit->set_id(key);
  }
  return it;
}

EICPIDParticle* EICPIDParticleContainer::findEICPIDParticle(EICPIDDefs::keytype key)
{
  EICPIDParticleContainer::ConstIterator it = m_particleMap.find(key);
  if (it != m_particleMap.end())
  {
    return it->second;
  }

  return nullptr;
}
