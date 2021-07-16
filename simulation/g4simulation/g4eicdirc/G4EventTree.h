#ifndef G4EVENTTREE_H
#define G4EVENTTREE_H

#define MAXHIT 500000

#include <Rtypes.h>

typedef struct
{
  // Event Level
  float energy;
  float theta;
  float phi;
  float px;
  float py;
  float pz;

  // Hit level
  int nhits;
  int detid[MAXHIT];
  int hitid[MAXHIT];
  int trkid[MAXHIT];
  float x0[MAXHIT];
  float y0[MAXHIT];
  float z0[MAXHIT];
  float x1[MAXHIT];
  float y1[MAXHIT];
  float z1[MAXHIT];
  float edep[MAXHIT];
  Int_t mcp_id[MAXHIT];
  Int_t pixel_id[MAXHIT];
  Double_t track_mom_bar[MAXHIT][3];
  Double_t track_hit_pos_bar[MAXHIT][3];
  Double_t lead_time[MAXHIT];
  Double_t wavelength[MAXHIT];
  Double_t hit_globalPos[MAXHIT][3];
  Double_t hit_localPos[MAXHIT][3];
  Double_t hit_digiPos[MAXHIT][3];
  Double_t hit_mom[MAXHIT][3];
  Double_t hit_pos[MAXHIT][3];
  Long64_t hit_pathId[MAXHIT];
  
} G4EventTree;

#endif
