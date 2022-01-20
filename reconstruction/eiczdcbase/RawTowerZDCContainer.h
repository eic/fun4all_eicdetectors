#ifndef EICZDCBASE_RAWTOWERZDCCONTAINER_H
#define EICZDCBASE_RAWTOWERZDCCONTAINER_H

#include "RawTowerZDCDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>

class RawTowerZDC;

class RawTowerZDCContainer : public PHObject
{
 public:
  typedef std::map<RawTowerZDCDefs::keytype, RawTowerZDC *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  RawTowerZDCContainer(RawTowerZDCDefs::CalorimeterId caloid = RawTowerZDCDefs::NONE)
    : _caloid(caloid)
  {
  }

  ~RawTowerZDCContainer() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream &os = std::cout) const override;

  void setCalorimeterID(RawTowerZDCDefs::CalorimeterId caloid) { _caloid = caloid; }
  RawTowerZDCDefs::CalorimeterId getCalorimeterID() { return _caloid; }

  ConstIterator AddTower(const unsigned int ieta, const unsigned int iphi, const unsigned int il, RawTowerZDC *twr);
  ConstIterator AddTower(RawTowerZDCDefs::keytype key, RawTowerZDC *twr);

  RawTowerZDC *getTower(RawTowerZDCDefs::keytype key);
  const RawTowerZDC *getTower(RawTowerZDCDefs::keytype key) const;

  RawTowerZDC *getTower(const unsigned int ieta, const unsigned int iphi, const unsigned int il );
  const RawTowerZDC *getTower(const unsigned int ieta, const unsigned int iphi, const unsigned int il) const;

  //! return all towers
  ConstRange getTowers(void) const;
  Range getTowers(void);

  unsigned int size() const { return _towers.size(); }
  void compress(const double emin);
  double getTotalEdep() const;

 protected:
  RawTowerZDCDefs::CalorimeterId _caloid;
  Map _towers;

  ClassDefOverride(RawTowerZDCContainer, 1)
};

#endif
