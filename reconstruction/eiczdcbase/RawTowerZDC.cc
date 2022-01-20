#include "RawTowerZDC.h"

#include <phool/phool.h>  // for PHOOL_VIRTUAL_WARN

#include <cstdlib>  // for exit

RawTowerZDC::CellMap DummyCellMap;
RawTowerZDC::ShowerMap DummyShowerMap;

RawTowerZDC::CellIterator RawTowerZDC::find_g4cell(CellKeyType /*id*/)
{
  return DummyCellMap.end();
}

RawTowerZDC::CellConstIterator RawTowerZDC::find_g4cell(CellKeyType /*id*/) const
{
  return DummyCellMap.end();
}

RawTowerZDC::CellConstRange RawTowerZDC::get_g4cells() const
{
  PHOOL_VIRTUAL_WARN("get_g4cells()");
  return CellConstRange(DummyCellMap.begin(), DummyCellMap.end());
}

RawTowerZDC::ShowerConstRange RawTowerZDC::get_g4showers() const
{
  PHOOL_VIRTUAL_WARN("get_g4showers()");
  return ShowerConstRange(DummyShowerMap.begin(), DummyShowerMap.end());
}

RawTowerZDC::ShowerIterator RawTowerZDC::find_g4shower(int /*id*/)
{
  return DummyShowerMap.end();
}

RawTowerZDC::ShowerConstIterator RawTowerZDC::find_g4shower(int /*id*/) const
{
  return DummyShowerMap.end();
}

const std::string RawTowerZDC::get_property_info(RawTowerZDC::PROPERTY prop_id)
{
  switch (prop_id)
  {
  case prop_scint_gammas:
    return "Scintillation photon count or energy";
  case prop_cerenkov_gammas:
    return "Cherenkov photon count or energy";

  default:
    std::cout << "RawTowerZDC::get_property_info - Fatal Error - unknown index " << prop_id << std::endl;
    exit(1);
  }
}
