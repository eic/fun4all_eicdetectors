#include "EICG4dRICHHit.h"

#include <phool/phool.h>

EICG4dRICHHit::EICG4dRICHHit(const PHG4Hit *g4hit) { CopyFrom(g4hit); };

void EICG4dRICHHit::Reset() 
{ // TODO: make sure this is filled out
  hitid = ULONG_LONG_MAX;
  trackid = INT_MIN;
  showerid = INT_MIN;
  edep = NAN;
  for (int i = 0; i < 2; i++) 
  {
    set_position(i, G4ThreeVector(NAN, NAN, NAN));
    set_t(i, NAN);
  };
};

int EICG4dRICHHit::get_detid() const 
{
  // a compile time check if the hit_idbits are within range (1-32)
  static_assert(PHG4HitDefs::hit_idbits <= sizeof(unsigned int) * 8,
                "hit_idbits < 32, fix in PHG4HitDefs.h");
  int detid = (hitid >> PHG4HitDefs::hit_idbits);
  return detid;
};

void EICG4dRICHHit::print() const 
{
  std::cout << "New EICG4dRICHHit  " << hitid << "  on track " << trackid << " EDep "
       << edep << std::endl;
  std::cout << "Location: X " << x[0] << "/" << x[1] << "  Y " << y[0] << "/" << y[1]
       << "  Z " << z[0] << "/" << z[1] << std::endl;
  std::cout << "Time        " << t[0] << "/" << t[1] << std::endl;
}
