#include "EICPIDParticle.h"

#include <TSystem.h>  // for gSystem

#include <cassert>
#include <cstdlib>
#include <type_traits>

using namespace std;

const std::map<EICPIDParticle::PROPERTY,
               std::pair<const std::string, EICPIDParticle::PROPERTY_TYPE> >
    EICPIDParticle::m_propertyInfo = {
        {Truth_PID, {"Truth PID", type_int}},
        {Truth_momentum, {"Truth Momentum", type_float}},
        {Truth_eta, {"Truth eta", type_float}},
        {CTTL_beta, {"Beta on CTTL", type_float}},
        {ETTL_beta, {"Beta on ETTL", type_float}},
        {FTTL_beta, {"Beta on FTTL", type_float}}  //
};

void EICPIDParticle::CopyFrom(const PHObject* phobj)
{
  const EICPIDParticle* src = dynamic_cast<const EICPIDParticle*>(phobj);
  assert(src);
  set_id(src->get_id());
  // This is a generic copy of ALL properties a hit has
  // do not add explicit copies, they will be added to
  // the new hits with their default value increasing memory use
  for (unsigned char ic = 0; ic < UCHAR_MAX; ic++)
  {
    PROPERTY prop_id = static_cast<EICPIDParticle::PROPERTY>(ic);
    if (src->has_property(prop_id))
    {
      set_property_nocheck(prop_id, src->get_property_nocheck(prop_id));
    }
  }
}

void EICPIDParticle::identify(ostream& os) const
{
  os << "Class " << this->ClassName() << endl;
  return;
}

void EICPIDParticle::Reset()
{
  cout << "Reset not implemented by daughter class" << endl;
  return;
}

std::pair<const std::string, EICPIDParticle::PROPERTY_TYPE>
EICPIDParticle::get_property_info(const PROPERTY prop_id)
{
  const auto iter = m_propertyInfo.find(prop_id);

  if (iter == m_propertyInfo.end())
  {
    cout << "EICPIDParticle::get_property_info - Fatal Error - unknown index " << prop_id << endl;
    gSystem->Exit(1);
    exit(1);
  }
  else
  {
    return iter->second;
  }
}

bool EICPIDParticle::check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type)
{
  pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
  if (property_info.second != prop_type)
  {
    return false;
  }
  return true;
}

string
EICPIDParticle::get_property_type(const PROPERTY_TYPE prop_type)
{
  switch (prop_type)
  {
  case type_int:
    return "int";
  case type_uint:
    return "unsigned int";
  case type_float:
    return "float";
  default:
    return "unkown";
  }
}
