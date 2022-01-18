#include "RawTowerZDCDeadMap.h"

#include <iostream>

const RawTowerZDCDeadMap::Map&
RawTowerZDCDeadMap::getDeadTowers(void) const
{
  static Map tmp_map;
  return tmp_map;
}

RawTowerZDCDeadMap::Map&
RawTowerZDCDeadMap::getDeadTowers(void)
{
  static Map tmp_map;
  return tmp_map;
}

void RawTowerZDCDeadMap::addDeadTower(const unsigned int /*ieta*/, const int unsigned /*iphi*/)
{
}

void RawTowerZDCDeadMap::addDeadTower(RawTowerZDCDefs::keytype /*key*/)
{
}

bool RawTowerZDCDeadMap::isDeadTower(RawTowerZDCDefs::keytype /*key*/)
{
  return false;
}

bool RawTowerZDCDeadMap::isDeadTower(const unsigned int /*ieta*/, const unsigned int /*iphi*/)
{
  return false;
}

int RawTowerZDCDeadMap::isValid() const
{
  return size() > 0;
}

void RawTowerZDCDeadMap::Reset()
{
}

void RawTowerZDCDeadMap::identify(std::ostream& os) const
{
  os << "RawTowerZDCDeadMap" << std::endl;
}
