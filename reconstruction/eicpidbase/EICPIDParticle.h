// TeLogLikelyhood emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICPID_EICPIDParticle_H
#define EICPID_EICPIDParticle_H

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include "EICPIDDefs.h"

class EICPIDParticle : public PHObject
{
 public:
  EICPIDParticle() {}
  ~EICPIDParticle() override {}

  void identify(std::ostream &os = std::cout) const override;
  void CopyFrom(const PHObject *phobj) override;
  void Reset() override;

  virtual EICPIDDefs::keytype get_id() const { return EICPIDDefs::INVALID_KEY; }
  virtual void set_id(const EICPIDDefs::keytype) { return; }

  virtual float get_SumLogLikelyhood(EICPIDDefs::PIDCandidate) const { return 0; }
  virtual float get_LogLikelyhood(EICPIDDefs::PIDCandidate, EICPIDDefs::PIDDetector) const { return m_minLogLikelihood; }
  virtual void set_LogLikelyhood(EICPIDDefs::PIDCandidate, EICPIDDefs::PIDDetector, float) {}

  enum PROPERTY_TYPE
  {  //
    type_int = 1,
    type_uint = 2,
    type_float = 3,
    type_unknown = -1
  };

  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to EICPIDParticle::m_propertyInfo
  enum PROPERTY
  {
    Truth_PID = 0,
    Truth_momentum ,
    Truth_eta,

    //
    CTTL_beta = 10,
    ETTL_beta,
    FTTL_beta,

    //! max limit in order to fit into 8 bit unsigned number
    prop_MAX_NUMBER = UCHAR_MAX
  };

  // property info definition in the EICPIDParticle.cc file
  const static std::map<PROPERTY,
                        std::pair<const std::string, EICPIDParticle::PROPERTY_TYPE> >
      m_propertyInfo;

  virtual bool has_property(const PROPERTY /*prop_id*/) const { return false; }
  virtual float get_property_float(const PROPERTY /*prop_id*/) const { return NAN; }
  virtual int get_property_int(const PROPERTY /*prop_id*/) const { return INT_MIN; }
  virtual unsigned int get_property_uint(const PROPERTY /*prop_id*/) const { return UINT_MAX; }
  virtual void set_property(const PROPERTY /*prop_id*/, const float /*value*/) { return; }
  virtual void set_property(const PROPERTY /*prop_id*/, const int /*value*/) { return; }
  virtual void set_property(const PROPERTY /*prop_id*/, const unsigned int /*value*/) { return; }
  static std::pair<const std::string, PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  static bool check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type);
  static std::string get_property_type(const PROPERTY_TYPE prop_type);

 protected:

  static constexpr float m_minLogLikelihood = -100;

  virtual unsigned int get_property_nocheck(const PROPERTY /*prop_id*/) const { return UINT_MAX; }
  virtual void set_property_nocheck(const PROPERTY /*prop_id*/, const unsigned int) { return; }
  ClassDefOverride(EICPIDParticle, 1)
};

#endif
