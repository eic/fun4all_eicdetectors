#include "RawTowerZDCGeomContainer.h"
#include "RawTowerZDCGeom.h"

#include <iostream>
#include <cassert>
#include <cstdlib>

using namespace std;

RawTowerZDCGeomContainer::RawTowerZDCGeomContainer(RawTowerZDCDefs::CalorimeterId caloid)
  : _caloid(caloid)
{
}

RawTowerZDCGeomContainer::~RawTowerZDCGeomContainer()
{
  Reset();  //make sure everything is deleted
}

void RawTowerZDCGeomContainer::identify(std::ostream& os) const
{
  os << "RawTowerZDCGeomContainer, number of tower geometries: " << size()
     << std::endl;
}

RawTowerZDCGeomContainer::ConstIterator
RawTowerZDCGeomContainer::add_tower_geometry(RawTowerZDCGeom* geo)
{
  assert(geo);

  if (RawTowerZDCDefs::decode_caloid(geo->get_id()) != get_calorimeter_id())
  {
    cout << "RawTowerZDCGeomContainer::add_tower_geometry - Fatal Error - "
            "attempting to add tower geometry with id = "
         << geo->get_id()
         << " with CaloID = " << RawTowerZDCDefs::decode_caloid(geo->get_id())
         << " to this container of CaloID = " << get_calorimeter_id() << ".";
    geo->identify(cout);
    exit(2);
  }

  Iterator it = _geoms.find(geo->get_id());
  if (it != _geoms.end())
  {
    cout
        << "RawTowerZDCGeomContainer::add_tower_geometry - WARNING - replace tower geometry for tower #"
        << geo->get_id() << ". This Old tower will be deleted: ";
    it->second->identify(cout);

    delete it->second;
    _geoms.erase(it);
  }

  _geoms[geo->get_id()] = geo;
  return _geoms.find(geo->get_id());
}

RawTowerZDCGeomContainer::ConstRange
RawTowerZDCGeomContainer::get_tower_geometries(void) const
{
  return make_pair<ConstIterator, ConstIterator>(_geoms.begin(), _geoms.end());
}

RawTowerZDCGeomContainer::Range
RawTowerZDCGeomContainer::get_tower_geometries(void)
{
  return make_pair<Iterator, Iterator>(_geoms.begin(), _geoms.end());
}

RawTowerZDCGeom*
RawTowerZDCGeomContainer::get_tower_geometry(RawTowerZDCDefs::keytype key)
{
  Iterator it = _geoms.find(key);
  if (it != _geoms.end())
  {
    return it->second;
  }
  return NULL;
}

int RawTowerZDCGeomContainer::isValid() const
{
  return (!_geoms.empty());
}

void RawTowerZDCGeomContainer::Reset()
{
  while (_geoms.begin() != _geoms.end())
  {
    delete _geoms.begin()->second;
    _geoms.erase(_geoms.begin());
  }
}
