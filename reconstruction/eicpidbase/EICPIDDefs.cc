#include "EICPIDDefs.h"

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>

namespace EICPIDDefs
{
PIDDetector getPIDDetector(const std::string& name)
{
  for (auto pair : PIDDetectorNameMap)
  {
    if (boost::iequals(pair.first, name))
      return pair.second;
  }
  return InvalidDetector;
}

const std::string& getPIDDetectorName(const PIDDetector det)
{
  static std::map<PIDDetector, std::string> reverse_map;

  if (reverse_map.size() == 0)
  {
    // build reverse map

    for (const auto& pair : PIDDetectorNameMap)
    {
      reverse_map[pair.second] = pair.first;
    }
  }

  if (reverse_map.find(det) == reverse_map.end())
  {
    std::cout << __PRETTY_FUNCTION__
              << " WARNING:  asking for det = " << det
              << ", which is not defined in EICPIDDefs::PIDDetectorNameMap" << std::endl;
  }

  return reverse_map[det];
}

}  // namespace EICPIDDefs
