#include "EICPIDDefs.h"

#include <boost/algorithm/string.hpp>
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

}  // namespace EICPIDDefs
