#include "RawTowerZDCContainer.h"

#include "RawTowerZDC.h"

#include <cstdlib>
#include <iostream>

using namespace std;

void RawTowerZDCContainer::compress(const double emin)
{
  if (emin <= 0)  // no need to loop through the map if we don't apply a cut
  {
    return;
  }
  Iterator itr = _towers.begin();
  Iterator last = _towers.end();
  for (; itr != last;)
  {
    RawTowerZDC *tower = (itr->second);
    if (tower->get_energy() < emin)
    {
      delete tower;
      _towers.erase(itr++);
    }
    else
    {
      ++itr;
    }
  }
}

RawTowerZDCContainer::ConstRange
RawTowerZDCContainer::getTowers(void) const
{
  return make_pair(_towers.begin(), _towers.end());
}

RawTowerZDCContainer::Range
RawTowerZDCContainer::getTowers(void)
{
  return make_pair(_towers.begin(), _towers.end());
}

RawTowerZDCContainer::ConstIterator
RawTowerZDCContainer::AddTower(const unsigned int ieta, const int unsigned iphi, const int unsigned il, RawTowerZDC *rawtower)
{
  RawTowerZDCDefs::keytype key = RawTowerZDCDefs::encode_towerid_zdc(_caloid, ieta, iphi, il);
  _towers[key] = rawtower;
  rawtower->set_id(key);  // force tower key to be synced to container key

  return _towers.find(key);
}

RawTowerZDCContainer::ConstIterator
RawTowerZDCContainer::AddTower(RawTowerZDCDefs::keytype key, RawTowerZDC *twr)
{
  if (RawTowerZDCDefs::decode_caloid(key) != _caloid)
  {
    cout << "RawTowerZDCContainer::AddTower - Error - adding tower to wrong container! Container CaloID = "
         << _caloid << ", requested CaloID = " << RawTowerZDCDefs::decode_caloid(key) << " based on key " << key << endl;
    exit(2);
  }

  _towers[key] = twr;
  twr->set_id(key);  // force tower key to be synced to container key

  return _towers.find(key);
}

RawTowerZDC *
RawTowerZDCContainer::getTower(RawTowerZDCDefs::keytype key)
{
  ConstIterator it = _towers.find(key);
  if (it != _towers.end())
  {
    return it->second;
  }
  return NULL;
}

const RawTowerZDC *
RawTowerZDCContainer::getTower(RawTowerZDCDefs::keytype key) const
{
  ConstIterator it = _towers.find(key);
  if (it != _towers.end())
  {
    return it->second;
  }
  return NULL;
}

RawTowerZDC *
RawTowerZDCContainer::getTower(const unsigned int ieta, const unsigned int iphi , const unsigned int il)
{
  RawTowerZDCDefs::keytype key = RawTowerZDCDefs::encode_towerid_zdc(_caloid, ieta, iphi, il);
  return getTower(key);
}

const RawTowerZDC *
RawTowerZDCContainer::getTower(const unsigned int ieta, const unsigned int iphi, const unsigned int il) const
{
  RawTowerZDCDefs::keytype key = RawTowerZDCDefs::encode_towerid_zdc(_caloid, ieta, iphi, il);
  return getTower(key);
}


int RawTowerZDCContainer::isValid() const
{
  return (!_towers.empty());
}

void RawTowerZDCContainer::Reset()
{
  while (_towers.begin() != _towers.end())
  {
    delete _towers.begin()->second;
    _towers.erase(_towers.begin());
  }
}

void RawTowerZDCContainer::identify(std::ostream &os) const
{
  os << "RawTowerZDCContainer, number of towers: " << size() << std::endl;
}

double
RawTowerZDCContainer::getTotalEdep() const
{
  double totalenergy = 0;
  ConstIterator iter;
  for (iter = _towers.begin(); iter != _towers.end(); ++iter)
  {
    totalenergy += iter->second->get_energy();
  }
  return totalenergy;
}
