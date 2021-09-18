// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICPID_EICPIDDefs_H
#define EICPID_EICPIDDefs_H

#include <climits>
#include <map>
#include <string>

namespace EICPIDDefs
{
// sync PID particle keys to the key of tracks
typedef unsigned int keytype;
static const keytype INVALID_KEY = UINT_MAX;

enum PIDDetector
{
  PIDAll = 0,
  mRICH = 1,
  DIRC = 2,
  dRICH = 3,
  GasRICH = 4,
  ETTL = 11,
  CTTL = 12,
  FTTL = 13,
  InvalidDetector = -1
};

const std::map<std::string, PIDDetector> PIDDetectorNameMap = {
    {"PIDAll", PIDAll},
    {"mRICH", mRICH},
    {"DIRC", DIRC},
    {"dRICH", dRICH},
    {"GasRICH", GasRICH},
    {"ETTL", ETTL},
    {"CTTL", CTTL},
    {"FTTL", FTTL}};

// consistent with PDG encoding
enum PIDCandidate
{
  ElectronCandiate = 11,
  MuonCandiate = 13,
  PionCandiate = 211,
  KaonCandiate = 321,
  ProtonCandiate = 2212,
  InvalidCandiate = 0
};

//! convert PID detector node names in to ID number for the container.
PIDDetector getPIDDetector(const std::string& name);

const std::string& getPIDDetectorName(const PIDDetector det);

}  // namespace EICPIDDefs

#endif
