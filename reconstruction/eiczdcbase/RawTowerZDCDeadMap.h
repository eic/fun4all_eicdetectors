#ifndef EICZDCBASE_RAWTOWERZDCDEADMAP_H
#define EICZDCBASE_RAWTOWERZDCDEADMAP_H

#include "RawTowerZDCDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <set>

class RawTowerZDCDeadMap : public PHObject
{
 public:
  typedef std::set<RawTowerZDCDefs::keytype> Map;

  ~RawTowerZDCDeadMap() override {}

  void Reset() override;
  int isValid() const override;

  void identify(std::ostream &os = std::cout) const override;

  virtual void setCalorimeterID(RawTowerZDCDefs::CalorimeterId /*caloid*/) {}
  virtual RawTowerZDCDefs::CalorimeterId getCalorimeterID() { return RawTowerZDCDefs::NONE; }
  virtual void addDeadTower(const unsigned int ieta, const unsigned int iphi);
  virtual void addDeadTower(RawTowerZDCDefs::keytype key);

  virtual bool isDeadTower(RawTowerZDCDefs::keytype key);
  virtual bool isDeadTower(const unsigned int ieta, const unsigned int iphi);
  //! return all towers
  virtual const Map &getDeadTowers(void) const;
  virtual Map &getDeadTowers(void);

  virtual unsigned int size() const { return 0; }

 protected:
  RawTowerZDCDeadMap(RawTowerZDCDefs::CalorimeterId /*caloid*/ = RawTowerZDCDefs::NONE)
  {
  }

 private:
  ClassDefOverride(RawTowerZDCDeadMap, 1)
};

#endif
