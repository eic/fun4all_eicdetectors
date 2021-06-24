#include "PHG4BarrelEcalDetector.h"
#include "PHG4BarrelEcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4DisplacedSolid.hh>
#include <Geant4/G4Tubs.hh>

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

using namespace std;

class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________________
PHG4BarrelEcalDetector::PHG4BarrelEcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4BarrelEcalDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_SupportActiveFlag(m_Params->get_int_param("supportactive"))
  , m_TowerLogicNamePrefix("bcalTower")
  , m_SuperDetector("NONE")
{
}
//_______________________________________________________________________
int PHG4BarrelEcalDetector::IsInBarrelEcal(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* mylogvol = volume->GetLogicalVolume();
  if (m_ActiveFlag)
  {
    if (m_ScintiLogicalVolSet.find(mylogvol) != m_ScintiLogicalVolSet.end())
    {
      return 1;
    }
  }

  if (m_AbsorberActiveFlag)
  {
    if (m_AbsorberLogicalVolSet.find(mylogvol) != m_AbsorberLogicalVolSet.end())
    {
      return -1;
    }
  }

  if (m_SupportActiveFlag)
  {
    if (m_SupportLogicalVolSet.find(mylogvol) != m_SupportLogicalVolSet.end())
    {
      return -2;
    }
  }
  return 0;
}

//_______________________________________________________________________
void PHG4BarrelEcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4BarrelEcalDetector: Begin Construction" << std::endl;
  }

  if (m_Params->get_string_param("mapping_file").empty())
  {
    std::cout << "ERROR in PHG4BarrelEcalDetector: No mapping file specified. Abort detector construction." << std::endl;
    std::cout << "Please run set_string_param(\"mapping_file\", std::string filename ) first." << std::endl;
    gSystem->Exit(1);
  }

  ParseParametersFromTable();

 const double radius = 85*cm;
 const double Length = 298.94*cm;
 const double max_radius = 138*cm;
 const double pos_x1 = 0*cm;
 const double pos_y1 = 0*cm;
 const double pos_z1 = 0*cm;

  G4Tubs *cylinder_solid = new G4Tubs("BCAL_SOLID",
                                       radius, max_radius,
                                       Length/ 2.0, 0, 2*M_PI);


  G4Material *cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);


  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid, cylinder_mat,
                                       "BCAL_SOLID", 0, 0, 0);

  m_DisplayAction->AddVolume(cylinder_logic, "BCalCylinder");

  //cylinder_physi = 

  std::string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  new G4PVPlacement(0, G4ThreeVector(pos_x1, pos_y1, pos_z1), cylinder_logic, name_envelope,
                                     logicWorld, false, 0, OverlapCheck());

  
  std::pair<G4LogicalVolume *, G4Transform3D> psec = Construct_AzimuthalSeg();
  G4LogicalVolume *sec_logic = psec.first;
  const G4Transform3D &sec_trans = psec.second;
  double Rot[32] = {0, 0.19635, 0.392699, 0.589049,0.785398, 0.981748, 1.1781, 1.37445, 1.5708, 1.76715, 1.9635, 2.15984, 2.35619, 2.55254, 2.74889, 2.94524, 3.14159, 3.33794, 3.53429, 3.73064, 3.92699, 4.12334, 4.31969, 4.51604, 4.71239, 4.90874, 5.10509, 5.30144, 5.49779, 5.69414, 5.89049, 6.08684};

  for (int i =0; i<32; i++)
  {
    const int sec = i;
    const double rot = Rot[i];

    G4Transform3D sec_place = G4RotateZ3D(rot) * sec_trans;

    stringstream name;
    
    name << GetName() << "_sec" << sec;

    new G4PVPlacement(sec_place, sec_logic,
                      "enclosure", cylinder_logic, false, sec,
                      OverlapCheck());
   }

  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  PlaceTower(cylinder_logic, singletower);

  return;
}


std::pair<G4LogicalVolume*, G4Transform3D>
PHG4BarrelEcalDetector::Construct_AzimuthalSeg()
{


  const double Max_radius =  138*cm;
  const double Thickness =  52.9998*cm;
  const double Radius = 85*cm;
  const int  phi_bin_in_sec = 4*cm;
  const double side_wall_thickness = 0.0762;
  const double sidewall_outer_torr = 0.15875;
  const double assembly_spacing = 0.01905;
  const double divider_width = 0;
  const double length = 298.94*cm;

  const int  Azimuthal_n_Sec = 256; 
  if (!(Azimuthal_n_Sec > 4))
  {
    cout << "azimuthal n sec <= 4: " << Azimuthal_n_Sec << endl;
    gSystem->Exit(1);
  }
  
  const G4double half_chord_backend =  Max_radius*tan(M_PI/Azimuthal_n_Sec) + fabs(Thickness*0.5);
  const G4double reduced_outer_radius = sqrt(pow(Max_radius, 2) - half_chord_backend * half_chord_backend);
  const G4double enclosure_depth = reduced_outer_radius - Radius;
  const G4double enclosure_center = 0.5 * (reduced_outer_radius + Radius); 
  const G4double enclosure_half_height_half_width = enclosure_center * tan(M_PI / Azimuthal_n_Sec);
  
  const G4double width_adj1 = tan(- M_PI / Azimuthal_n_Sec) * enclosure_depth * 0.5;
  const G4double width_adj2 = tan(+ M_PI / Azimuthal_n_Sec) * enclosure_depth * 0.5;

  const G4double center_adj = (width_adj1 + width_adj2) * 0.5;
  const G4double center_tilt_angle = atan2(center_adj, enclosure_depth * 0.5);
  
  // enclosure walls
  const G4double edge1_tilt_angle = atan2(width_adj1, enclosure_depth * 0.5);
  const G4double edge2_tilt_angle = atan2(width_adj2, enclosure_depth * 0.5);
  
  // projective center
  const G4double half_projection_ratio = 0.5 * (-width_adj1 + width_adj2) / enclosure_half_height_half_width;
  const G4double projection_center_y = enclosure_center - ((enclosure_depth * 0.5) / half_projection_ratio);
  const G4double projection_center_x = center_adj / half_projection_ratio;

  // blocks azimuthal segmentation 
  assert(phi_bin_in_sec >= 1);
  const G4double block_azimuth_angle = (edge2_tilt_angle - edge1_tilt_angle) / phi_bin_in_sec;
  assert(block_azimuth_angle > 0);

  if (!(fabs(block_azimuth_angle - M_PI * 2 / Azimuthal_n_Sec / phi_bin_in_sec) < M_PI * numeric_limits<G4double>::epsilon()))
  {
    cout << "angle/nsec out of range: " << M_PI * numeric_limits<G4double>::epsilon() << endl;
    gSystem->Exit(1);
  }
  const G4double block_edge1_half_width = enclosure_half_height_half_width - (side_wall_thickness+sidewall_outer_torr + 2.0 * assembly_spacing) / cos(edge1_tilt_angle);
  const G4double block_edge2_half_width = enclosure_half_height_half_width - (side_wall_thickness+sidewall_outer_torr + 2.0 * assembly_spacing ) / cos(edge2_tilt_angle);

  G4double block_width_ratio = 0;
  for (int s = 0; s < phi_bin_in_sec; ++s)
  {
    block_width_ratio += 1 / cos(block_azimuth_angle * (0.5 + s) + edge1_tilt_angle);

  }
  const G4double block_half_height_width =(block_edge1_half_width + block_edge2_half_width) / block_width_ratio;

  assert(block_half_height_width > 0);

 // write out the azimuthal block geometry
  // block azimuth geometry records 
  
  struct block_azimuth_geom
  {
    G4double angle;
    G4double projection_center_y;
    G4double projection_center_x;
    G4double projection_length;
  };
  
  vector<block_azimuth_geom> block_azimuth_geoms(phi_bin_in_sec,
                                                 block_azimuth_geom{
                                                     numeric_limits<double>::signaling_NaN(),
                                                     numeric_limits<double>::signaling_NaN(),
                                                     numeric_limits<double>::signaling_NaN(),
                                                     numeric_limits<double>::signaling_NaN()});  // [phi-bin in sector] -> azimuth geometry
  
  G4double block_x_edge1 = block_edge1_half_width;


  for (int s = 0; s < phi_bin_in_sec; ++s)
  {
    block_azimuth_geom& geom = block_azimuth_geoms[s];

    geom.angle = block_azimuth_angle * (0.5 + s) + edge1_tilt_angle;

    const G4double block_x_size = block_half_height_width / cos(geom.angle);
    assert(block_x_size > 0);
    const G4double x_center = block_x_edge1 - 0.5 * block_x_size;
    // projection center per block
    geom.projection_length = block_half_height_width / 2. / tan(block_azimuth_angle / 2.);
    assert(geom.projection_length > 0);
    geom.projection_center_y = enclosure_center - geom.projection_length * cos(geom.angle);
    geom.projection_center_x = x_center + geom.projection_length * sin(geom.angle);
    // next step
    block_x_edge1 -= block_x_size;

  }

  //write out the azimuthal block divider's geometry
  struct block_divider_azimuth_geom
  {
    G4double angle;  //! rotation angle
    G4double projection_center_y;
    G4double projection_center_x;
    G4double thickness;            // thickness in the approximate azimuth direction
    G4double radial_displacement;  //! displacement along the width direction, which is the radial direction if tilt = 0
    G4double width;                //! wdith along the approximate radial direction
  };
  assert(phi_bin_in_sec >= 1);
  vector<block_divider_azimuth_geom> divider_azimuth_geoms(phi_bin_in_sec - 1,
                                                           block_divider_azimuth_geom{
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN()});

  if (side_wall_thickness> 0 )
  {
    for (int s = 0; s < phi_bin_in_sec - 1; ++s)
    {
      block_divider_azimuth_geom& geom = divider_azimuth_geoms[s];

      geom.angle = 0.5 * (block_azimuth_geoms[s].angle + block_azimuth_geoms[s + 1].angle);
    //  geom.projection_center_y = 0.5 * (block_azimuth_geoms[s].projection_center_y + block_azimuth_geoms[s + 1].projection_center_y);

      geom.projection_center_y = enclosure_center;
      geom.projection_center_x = 0.5 * (block_azimuth_geoms[s].projection_center_x + block_azimuth_geoms[s + 1].projection_center_x);
      geom.radial_displacement = 0.5 * (block_azimuth_geoms[s].projection_length + block_azimuth_geoms[s + 1].projection_length);

      geom.thickness = 2.0 *assembly_spacing*cos(block_azimuth_angle / 2.) - 2 * um;
      geom.width = divider_width;
    }
  }

  if (fabs(block_x_edge1 - (-block_edge2_half_width)) > assembly_spacing )
  {
    cout << "PHG4FullProjBarrelCaloDetector2::Construct_AzimuthalSeg - ERROR - " << endl
         << "\t block_x_edge1 = " << block_x_edge1 << endl
         << "\t block_edge2_half_width = " << block_edge2_half_width << endl
         << "\t fabs(block_x_edge1 - (-block_edge2_half_width)) = " << fabs(block_x_edge1 - (-block_edge2_half_width)) << endl
         << "\t assembly_spacing() " << assembly_spacing  << endl;
    cout << "closure check failed: " << fabs(block_x_edge1 - (-block_edge2_half_width)) << endl;
    gSystem->Exit(1);
  }

  if (Verbosity())
  {
    cout << "PHG4FullProjBarrelCaloDetector2::Construct_AzimuthalSeg - " << endl
         << "\t edge1_tilt_angle = " << edge1_tilt_angle << endl
         << "\t edge2_tilt_angle = " << edge2_tilt_angle << endl
         << "\t projection_center_y = " << projection_center_y << endl
         << "\t projection_center_x = " << projection_center_x << endl
         << "\t block_azimuth_angle = " << block_azimuth_angle << endl
         << "\t block_edge1_half_width = " << block_edge1_half_width << endl
         << "\t block_edge2_half_width = " << block_edge2_half_width << endl
         << "\t block_width_ratio = " << block_width_ratio << endl
         << "\t block_half_height_width = " << block_half_height_width << endl;

    for (int s = 0; s < phi_bin_in_sec; ++s)
    {
      cout << "\t block[" << s << "].angle = " << block_azimuth_geoms[s].angle << endl;
      cout << "\t block[" << s << "].projection_center_y = " << block_azimuth_geoms[s].projection_center_y << endl;
      cout << "\t block[" << s << "].projection_center_x = " << block_azimuth_geoms[s].projection_center_x << endl;
    }
    for (int s = 0; s < phi_bin_in_sec - 1; ++s)
    {
      cout << "\t divider[" << s << "].angle = " << divider_azimuth_geoms[s].angle << endl;
      cout << "\t divider[" << s << "].projection_center_x = " << divider_azimuth_geoms[s].projection_center_x << endl;
      cout << "\t divider[" << s << "].projection_center_y = " << divider_azimuth_geoms[s].projection_center_y << endl;
      cout << "\t divider[" << s << "].radial_displacement = " << divider_azimuth_geoms[s].radial_displacement << endl;
      cout << "\t divider[" << s << "].thickness = " << divider_azimuth_geoms[s].thickness << endl;
      cout << "\t divider[" << s << "].width = " << divider_azimuth_geoms[s].width << endl;
    }
  }

  double halfpi = M_PI/2;

  assert(enclosure_depth > 10 * cm);

  G4VSolid* sec_solid = new G4Trap(
      G4String(GetName() + string("_sec_trap")),
      enclosure_depth * 0.5,                                                                              // G4double pDz,
      center_tilt_angle, halfpi,                                                                          // G4double pTheta, G4double pPhi
      //inner_half_width, length / 2.0, length  / 2.0,   // G4double pDy1, G4double pDx1, G4double pDx2,
      0.01, length / 2.0, length  / 2.0,   // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                                  // G4double pAlp1,
      //outter_half_width, length / 2.0, length / 2.0,  // G4double pDy2, G4double pDx3, G4double pDx4,
      0.01,  length / 2.0, length / 2.0,
      0                                                                                                   // G4double pAlp2 //
  );

 
  G4Transform3D sec_solid_transform =
      G4TranslateY3D(enclosure_center) * G4RotateY3D(halfpi) * G4RotateX3D(-halfpi);

  G4VSolid* sec_solid_place = new G4DisplacedSolid(
      G4String(GetName() + string("_sec")), sec_solid, sec_solid_transform);

  G4Material* cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);

  G4LogicalVolume* sec_logic = new G4LogicalVolume(sec_solid_place, cylinder_mat,
                                                   G4String(G4String(GetName() + string("_sec"))), 0, 0, nullptr);

  m_DisplayAction->AddVolume(sec_logic, "Sector");

  m_SupportLogicalVolSet.insert(sec_logic);
  return make_pair(sec_logic, G4Transform3D::Identity);

}


G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructTower()
{

  G4Trap* block_solid = new G4Trap(
                  "solid_tower",
      22.75*cm,                                                 // G4double pDz,
      0,  0,                                                    // G4double pTheta, G4double pPhi,
      2.0*cm, 2.0*cm,2.0*cm,                                    // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                       // G4double pAlp1,
      2.5*cm, 2.5*cm,2.5*cm,                                    // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                         // G4double pAlp2 //
  );


  G4Material* cylinder_mat = G4Material::GetMaterial("sciglass");
  assert(cylinder_mat);

  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, cylinder_mat,
                                                     "solid_tower", 0, 0,
                                                     nullptr);
  m_ScintiLogicalVolSet.insert(block_logic);
  m_DisplayAction->AddVolume(block_logic, "Block");

  return block_logic;
}

int PHG4BarrelEcalDetector::PlaceTower(G4LogicalVolume* sec, G4LogicalVolume* singletower)
{
  /* Loop over all tower positions in vector and place tower */
  for (std::map<std::string, towerposition>::iterator iterator = m_TowerPostionMap.begin(); iterator != m_TowerPostionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4BarrelEcalDetector: Place tower " << iterator->first
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k<< std::endl;
                //<< " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << std::endl;
    }
     int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;

    G4Transform3D block_trans =
        G4TranslateX3D(iterator->second.centerx) *
        G4TranslateY3D(iterator->second.centery) *
        G4RotateZ3D(iterator->second.azangle) *
        G4TranslateX3D(iterator->second.towercenterx) *
        G4TranslateY3D(iterator->second.towercentery) *
        G4TranslateZ3D(iterator->second.towercenterz) *
        G4RotateX3D(iterator->second.rotx);

    G4Transform3D sec_place = G4RotateZ3D(iterator->second.rot) * block_trans;   

    
    
    new G4PVPlacement(sec_place,
                      singletower,
                      iterator->first,
                      sec,
                      0, copyno, OverlapCheck());
    
  }

  return 0;
}


int PHG4BarrelEcalDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  std::ifstream istream_mapping;
  istream_mapping.open(m_Params->get_string_param("mapping_file"));
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in PHG4BarrelEcalDetector: Failed to open mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
    gSystem->Exit(1);
  }

  /* loop over lines in file */
  std::string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {

    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != std::string::npos)
    {
      if (Verbosity() > 0)
      {
        std::cout << "PHG4BarrelEcalDetector: SKIPPING line in mapping file: " << line_mapping << std::endl;
      }
      continue;
    }

    std::istringstream iss(line_mapping);

      unsigned idphi_j, ideta_k;
      G4double cx, cy, aa;
      G4double tcx, tcy, tcz;
      G4double rot_x, rott;
      std::string dummys;
 
       if (!(iss >> dummys >> idphi_j >> ideta_k >> cx >> cy >> aa >> tcx >> tcy >> tcz >>  rot_x >>  rott))
      {
        std::cout << "ERROR in PHG4BarrelEcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
        gSystem->Exit(1);
      }

    /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      std::ostringstream towername;
      towername.str("");
      towername << m_TowerLogicNamePrefix << "_j_" << idphi_j << "_k_" << ideta_k;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.centerx = cx;
      tower_new.centery = cy;
      tower_new.azangle = aa;
      tower_new.towercenterx = tcx;
      tower_new.towercentery = tcy;
      tower_new.towercenterz = tcz;
      tower_new.rotx = rot_x;
      tower_new.rot = rott;
      tower_new.idx_j = idphi_j;
      tower_new.idx_k = ideta_k;
      m_TowerPostionMap.insert(make_pair(towername.str(), tower_new));

      
    }

  return 0;
}