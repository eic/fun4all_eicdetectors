#include "EICPIDParticlev1.h"

#include <phool/phool.h>

#include <climits>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include "EICPIDDefs.h"

using namespace std;

EICPIDParticlev1::EICPIDParticlev1(const EICPIDParticle* g4hit)
{
  CopyFrom(g4hit);
}

void EICPIDParticlev1::Reset()
{
  m_id = -1;
  prop_map.clear();
}

bool EICPIDParticlev1::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i != prop_map.end();
}

float EICPIDParticlev1::get_property_float(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_float))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_float) << endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).fdata;

  return NAN;
}

int EICPIDParticlev1::get_property_int(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_int))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_int) << endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).idata;

  return INT_MIN;
}

unsigned int
EICPIDParticlev1::get_property_uint(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_uint))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_uint) << endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).uidata;

  return UINT_MAX;
}

void EICPIDParticlev1::set_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id, type_float))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_float) << endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

void EICPIDParticlev1::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id, type_int))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_int) << endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

void EICPIDParticlev1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id, type_uint))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_uint) << endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

unsigned int
EICPIDParticlev1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
  {
    return iter->second;
  }
  return UINT_MAX;
}

void EICPIDParticlev1::identify(ostream& os) const
{
  os << "Class " << this->ClassName();
  os << " id: " << m_id << endl;
  for (prop_map_t::const_iterator i = prop_map.begin(); i != prop_map.end(); ++i)
  {
    PROPERTY prop_id = static_cast<PROPERTY>(i->first);
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    os << "\t" << prop_id << ":\t" << property_info.first << " = \t";
    switch (property_info.second)
    {
    case type_int:
      os << get_property_int(prop_id);
      break;
    case type_uint:
      os << get_property_uint(prop_id);
      break;
    case type_float:
      os << get_property_float(prop_id);
      break;
    default:
      os << " unknown type ";
    }
    os << endl;
  }
}
