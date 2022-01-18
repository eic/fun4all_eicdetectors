#ifndef EICZDCBASE_RAWTOWERZDCGEOMCONTAINER_H
#define EICZDCBASE_RAWTOWERZDCGEOMCONTAINER_H

#include "RawTowerZDCDefs.h"

#include <phool/PHObject.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <type_traits>
#include <utility>


class RawTowerZDCGeom;

/*! \class RawTowerZDCGeomContainer
    \brief  Generic tower geometry class, store each tower's geometry individually
*/
class RawTowerZDCGeomContainer : public PHObject
{
 public:
  typedef std::map<RawTowerZDCDefs::keytype, RawTowerZDCGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  //! default constructor for ROOT IO
  RawTowerZDCGeomContainer(RawTowerZDCDefs::CalorimeterId caloid = RawTowerZDCDefs::NONE);
  ~RawTowerZDCGeomContainer() override;

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream &os = std::cout) const override;


  void set_calorimeter_id(RawTowerZDCDefs::CalorimeterId caloid) { _caloid = caloid; }
  RawTowerZDCDefs::CalorimeterId get_calorimeter_id() { return _caloid; }

  //! go through all towers
  ConstIterator add_tower_geometry(RawTowerZDCGeom *geo) ;
  RawTowerZDCGeom *get_tower_geometry(RawTowerZDCDefs::keytype key);

  //! return all tower geometries

  ConstRange get_tower_geometries(void) const;
  Range get_tower_geometries(void);

  unsigned int size() const { return _geoms.size(); }

 protected:
  RawTowerZDCDefs::CalorimeterId _caloid;
  Map _geoms;
  
  ClassDefOverride(RawTowerZDCGeomContainer, 1)
};

#endif
