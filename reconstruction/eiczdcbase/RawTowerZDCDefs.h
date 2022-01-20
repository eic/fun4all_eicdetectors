#ifndef EICZDCBASE_RAWTOWERZDCDEFS_H
#define EICZDCBASE_RAWTOWERZDCDEFS_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <bitset>

/*! Namespace with functions to encode / decode CaloTowerID. The highest 8 bits of the tower ID encode a unique ID
 * for the calorimeter the tower is in. The lower 24 bits uniquely identify the tower within a calorimeter.
 *
 */
namespace RawTowerZDCDefs
{
  /*! Define data type of unique tower ID, i.e. for CaloTowerID
   */
  typedef unsigned int keytype;

  /*! Bit ranges for encoding calorimeter ID and tower indices in combined tower ID
   */
  static unsigned int calo_idbits = 4;
  static unsigned int tower_idbits = sizeof(keytype) * 8 - calo_idbits;
  //  static unsigned int bitsIndex1 = 10;   //max 0x3FF (1023)
  static unsigned int bitsIndex2 = 10; // max 0x3FF (1023)
  static unsigned int bitsIndex3 = 8;  // max 0xFF (255)
  static unsigned int maxbitsCaloId = 0xF;
  static unsigned int maxbitsIndex1 = 0x3FF;
  static unsigned int maxbitsIndex2 = 0x3FF;
  static unsigned int maxbitsIndex3 = 0xFF;

  /*! Enum with all available calorimeter IDs. This enum can be extended up to 16 entries.
   * If adding new CalorimeterIDs, please also add them to the decode_caloname function below.
   */
  enum CalorimeterId
  {
    NONE = 0,
    ZDC_Crystal = 1,
    ZDC_SiPixel = 2,
    ZDC_SiPad =3,
    ZDC_Sci = 4,
  };

  /*! Extract calorimeter ID from CaloTowerID
   */
  inline CalorimeterId
  decode_caloid(const unsigned int calo_tower_id)
  {
    return static_cast<CalorimeterId>((calo_tower_id >> RawTowerZDCDefs::tower_idbits) & maxbitsCaloId);
  }

   /*! Extract tower index 1 of calorimeter tower from CaloTowerID with 3 indices for ZDC
   */
  inline unsigned int
  decode_index1zdc(const unsigned int calo_tower_id)
  {

    return (calo_tower_id >> (bitsIndex2+bitsIndex3)) & maxbitsIndex1;
  }

   /*! Extract tower index 2 of calorimeter tower from CaloTowerID with 3 indices for ZDC
   */
  inline unsigned int
  decode_index2zdc(const unsigned int calo_tower_id)
  {

    return (calo_tower_id >> (bitsIndex3)) & maxbitsIndex2;
  }

  /*! Extract tower index 3 of calorimeter tower from CaloTowerID with 3 indices for ZDC
   */
  inline unsigned int
  decode_index3zdc(const unsigned int calo_tower_id)
  {

    return calo_tower_id & maxbitsIndex3;
  }


   /*! Returns CaloTowerID for given calorimeter ID, tower index 1, tower index 2 and tower index 3
   */
  inline RawTowerZDCDefs::keytype
  encode_towerid_zdc(const CalorimeterId calo_id, 
		     const unsigned int tower_index_1,
		     const unsigned int tower_index_2, 
		     const unsigned int tower_index_3)
  {
    RawTowerZDCDefs::keytype calo_tower_id = 0;
    
    if (tower_index_1 < maxbitsIndex1 && tower_index_2 < maxbitsIndex2 &&  tower_index_3 < maxbitsIndex3)
    {
      calo_tower_id = (calo_id << RawTowerZDCDefs::tower_idbits) + (tower_index_1 << (bitsIndex2+bitsIndex3)) + (tower_index_2 << bitsIndex3) + tower_index_3;
    }
    else
    {
      std::cout << "too large index1 and/or index2; index1: "
                << tower_index_1 << " (max val " << maxbitsIndex1 << ")"
                << ", index2: "
                << tower_index_2 << " (max val " << maxbitsIndex2 << ")"
                << ", index3: "
                << tower_index_3 << " (max val " << maxbitsIndex3 << ")"<< std::endl;
      exit(1);
    }

    return calo_tower_id;
  }

  /*! Convert calorimeter ID to name string
   */
  inline std::string
  convert_caloid_to_name(const RawTowerZDCDefs::CalorimeterId calo_id)
  {
    switch (calo_id)
    {
    case NONE:
      return "NONE";
      break;

    case ZDC_Crystal:
      return "ZDC_Crtstal";
      break;
  
    case ZDC_SiPixel:
      return "ZDC_SiPixel";
      break;
  
    case ZDC_SiPad:
      return "ZDC_SiPad";
      break;
  
    case ZDC_Sci:
      return "ZDC_Sci";
      break;
  
    default:
      std::cout
          << "Invalid calorimeter ID passed to RawTowerZDCDefs::convert_caloid_to_name"
          << std::endl;
      exit(1);
    }
  }

  /*! Convert name string to calorimeter ID
   */
  inline RawTowerZDCDefs::CalorimeterId
  convert_name_to_caloid(const std::string &caloname)
  {
    if (caloname == "NONE")
      return NONE;

    else if (caloname == "ZDC_Crystal")
      return ZDC_Crystal;

    else if (caloname == "ZDC_SiPixel")
      return ZDC_SiPixel;

    else if (caloname == "ZDC_SiPad")
      return ZDC_SiPad;

    else if (caloname == "ZDC_Sci")
      return ZDC_Sci;

    else
    {
      std::cout << "Invalid calorimeter name " << caloname
                << " passed to RawTowerZDCDefs::convert_name_to_caloid" << std::endl;
      exit(1);
    }
  }

}  // end namespace RawTowerZDCDefs

#endif
