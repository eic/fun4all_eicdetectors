#include "RawTowerZDCv1.h"

#include <cmath>
#include <iostream>

using namespace std;

RawTowerZDCv1::RawTowerZDCv1(const RawTowerZDC& tower)
{
  towerid = (tower.get_id());
  energy = (tower.get_energy());
  time = (tower.get_time());

  CellConstRange cell_range = tower.get_g4cells();

  for (CellConstIterator cell_iter = cell_range.first;
       cell_iter != cell_range.second; ++cell_iter)
  {
    add_ecell(cell_iter->first, cell_iter->second);
  }

  ShowerConstRange shower_range = tower.get_g4showers();

  for (ShowerConstIterator shower_iter = shower_range.first;
       shower_iter != shower_range.second; ++shower_iter)
  {
    add_eshower(shower_iter->first, shower_iter->second);
  }
}

RawTowerZDCv1::RawTowerZDCv1(RawTowerZDCDefs::keytype id)
  : towerid(id)
{
}

void RawTowerZDCv1::Reset()
{
  energy = 0;
  time = NAN;
  ecells.clear();
  eshowers.clear();
}

int RawTowerZDCv1::isValid() const
{
  return get_energy() != 0;
}

void RawTowerZDCv1::identify(std::ostream& os) const
{
  
  os << "RawTowerZDCv1: etabin: " << get_bineta() << ", phibin: " << get_binphi() << ", l-bin: " << get_binl()
     << " energy=" << get_energy() << std::endl;
  return;
    
}

void RawTowerZDCv1::add_ecell(const CellKeyType g4cellid,
                           const float ecell)
{
  if (ecells.find(g4cellid) == ecells.end())
  {
    ecells[g4cellid] = ecell;
  }
  else
  {
    ecells[g4cellid] += ecell;
  }
}

void RawTowerZDCv1::add_eshower(const int g4showerid, const float eshower)
{
  if (eshowers.find(g4showerid) == eshowers.end())
  {
    eshowers[g4showerid] = eshower;
  }
  else
  {
    eshowers[g4showerid] += eshower;
  }
}

int RawTowerZDCv1::get_bineta() const
{
    return RawTowerZDCDefs::decode_index1zdc(towerid); 
}

int RawTowerZDCv1::get_binphi() const
{
  return RawTowerZDCDefs::decode_index2zdc(towerid);
}

int RawTowerZDCv1::get_binl() const
{
  return RawTowerZDCDefs::decode_index3zdc(towerid);
  
}

