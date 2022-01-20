#ifndef EICZDCBASE_RAWTOWERV1_H
#define EICZDCBASE_RAWTOWERV1_H

#include "RawTowerZDC.h"

#include "RawTowerZDCDefs.h"

#include <cstddef>
#include <iostream>
#include <map>
#include <utility>

class RawTowerZDCv1 : public RawTowerZDC
{
 public:
  RawTowerZDCv1(){}
  RawTowerZDCv1(const RawTowerZDC& tower);
  RawTowerZDCv1(RawTowerZDCDefs::keytype id);
  ~RawTowerZDCv1() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  void set_id(RawTowerZDCDefs::keytype id) override { towerid = id; }
  RawTowerZDCDefs::keytype get_id() const override { return towerid; }
  int get_bineta() const override;
  int get_binphi() const override;
  int get_binl() const  override;

  double get_energy() const override { return energy; }
  void set_energy(const double e) override { energy = e; }
  float get_time() const override { return time; }
  void set_time(const float t) override { time = t; }

  //---cell access--------------------------------------------------------------

  bool empty_g4cells() const override { return ecells.empty(); }
  size_t size_g4cells() const override { return ecells.size(); }
  RawTowerZDC::CellConstRange get_g4cells() const override
  {
    return make_pair(ecells.begin(), ecells.end());
  }
  RawTowerZDC::CellIterator find_g4cell(CellKeyType id) override { return ecells.find(id); }
  RawTowerZDC::CellConstIterator find_g4cell(CellKeyType id) const override { return ecells.find(id); }
  void add_ecell(const CellKeyType g4cellid,
                 const float ecell) override;
  void clear_g4cells() override { ecells.clear(); }

  //---shower access------------------------------------------------------------

  bool empty_g4showers() const override { return eshowers.empty(); }
  size_t size_g4showers() const override { return eshowers.size(); }
  RawTowerZDC::ShowerConstRange get_g4showers() const override
  {
    return make_pair(eshowers.begin(), eshowers.end());
  }
  RawTowerZDC::ShowerIterator find_g4shower(int id) override { return eshowers.find(id); }
  RawTowerZDC::ShowerConstIterator find_g4shower(int id) const override { return eshowers.find(id); }
  void add_eshower(const int g4showerid, const float eshower) override;
  void clear_g4showers() override { eshowers.clear(); }

 protected:
  RawTowerZDCDefs::keytype towerid = ~0;

  //! energy assigned to the tower. Depending on stage of process and DST node
  //! name, it could be energy deposition, light yield or calibrated energies
  double energy = 0.;
  //! Time stamp assigned to the tower. Depending on the tower maker, it could
  //! be rise time or peak time.
  float time = NAN;

  CellMap ecells;      //< default truth storage
  ShowerMap eshowers;  //< alternate truth storage for smaller filesizes

  ClassDefOverride(RawTowerZDCv1, 1)
};

#endif
