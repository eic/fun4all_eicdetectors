#include "EICPIDDefs.h"

#include <tr1/functional>

namespace EICPIDDefs
{

  int get_volume_id(const std::string & nodename)
  {
    return std::tr1::hash<std::string>()(nodename);
  }

}


