// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICPID_EICPIDParticleV1_H
#define EICPID_EICPIDParticleV1_H

#include <climits>  // for INT_MIN, ULONG_LONG_MAX
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>

#include "EICPIDDefs.h"
#include "EICPIDParticle.h"

class EICPIDParticlev1 : public EICPIDParticle
{
 public:
  EICPIDParticlev1() = default;
  explicit EICPIDParticlev1(const EICPIDParticle* g4hit);
  ~EICPIDParticlev1() override = default;
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;

  EICPIDDefs::keytype get_id() const override { return m_id; }
  void set_id(const EICPIDDefs::keytype i) override { m_id = i; }


  bool has_property(const PROPERTY prop_id) const override;
  float get_property_float(const PROPERTY prop_id) const override;
  int get_property_int(const PROPERTY prop_id) const override;
  unsigned int get_property_uint(const PROPERTY prop_id) const override;
  void set_property(const PROPERTY prop_id, const float value) override;
  void set_property(const PROPERTY prop_id, const int value) override;
  void set_property(const PROPERTY prop_id, const unsigned int value) override;

 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const override;
  void set_property_nocheck(const PROPERTY prop_id, const unsigned int ui) override { prop_map[prop_id] = ui; }

  EICPIDDefs::keytype m_id = -1;

  //! storage types for additional property
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;

  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property {
    float fdata;
    int32_t idata;
    uint32_t uidata;

    u_property(int32_t in)
      : idata(in)
    {
    }
    u_property(uint32_t in)
      : uidata(in)
    {
    }
    u_property(float in)
      : fdata(in)
    {
    }
    u_property()
      : uidata(0)
    {
    }

    prop_storage_t get_storage() const { return uidata; }
  };

  //! container for additional property
  prop_map_t prop_map;

  ClassDefOverride(EICPIDParticlev1, 2)
};

#endif
