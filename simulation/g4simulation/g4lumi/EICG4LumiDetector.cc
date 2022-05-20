//____________________________________________________________________________..
//
// This is a working template for the G4 Construct() method which needs to be implemented
// We wedge a method between the G4 Construct() to enable volume hierarchies on the macro
// so here it is called ConstructMe() but there is no functional difference
// Currently this installs a simple G4Box solid, creates a logical volume from it
// and places it. Put your own detector in place (just make sure all active volumes
// get inserted into the m_PhysicalVolumesSet)
// 
// Rather than using hardcoded values you should consider using the parameter class
// Parameter names and defaults are set in EICG4LumiSubsystem::SetDefaultParameters()
// Only parameters defined there can be used (also to override in the macro)
// to avoids typos.
// IMPORTANT: parameters have no inherent units, there is a convention (cm/deg)
// but in any case you need to multiply them here with the correct CLHEP/G4 unit 
// 
// The place where you put your own detector is marked with
// //begin implement your own here://
// //end implement your own here://
// Do not forget to include the G4 includes for your volumes
//____________________________________________________________________________..

#include "EICG4LumiDetector.h"


#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Color.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <TSystem.h>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <phool/recoConsts.h> //For rc WorldMaterial

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <utility>

class G4VSolid;
class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
EICG4LumiDetector::EICG4LumiDetector(PHG4Subsystem *subsys,
                                         PHCompositeNode *Node,
                                         PHParameters *parameters,
                                         const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  ,m_Layer(lyr)

{
}

int EICG4LumiDetector::IsInDetector(G4VPhysicalVolume *volume) const
{

  if( m_ActivePhysicalVolumesMap.find( volume ) != m_ActivePhysicalVolumesMap.end() )
  {
    return 1;
  }

  return 0;

}

//_______________________________________________________________
int EICG4LumiDetector::IsInVirtualDetector(G4VPhysicalVolume *volume) const
{

  if( m_VirtualPhysicalVolumesMap.find( volume ) != m_VirtualPhysicalVolumesMap.end() )
  {
    return 1;
  }

  return 0;

}

//_______________________________________________________________
int EICG4LumiDetector::GetDetId(G4VPhysicalVolume *volume) const
{
  
  if (IsInDetector(volume))
  {
    return 1;
  }

  return -1;
}

//_______________________________________________________________
void EICG4LumiDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
 //begin implement your own here://
 // Do not forget to multiply the parameters with their respective CLHEP/G4 unit !
  if (Verbosity() > 0)
  {
    std::cout << "EICG4LumiDetector: Begin Construction" << std::endl;
  }

  SetParametersFromFile();

  double enclosureCenter = m_Params->get_double_param( "FBenclosure_center" ) * cm;
  double LumiSpec_Z = m_Params->get_double_param( "LumiSpec_Z" ) * cm - enclosureCenter;
  double LumiSpec_XY = m_Params->get_double_param( "LumiSpec_XY" ) * cm;
  double LumiPhotonCAL_Z = m_Params->get_double_param( "LumiPhotonCAL_Z" ) * cm - enclosureCenter;
  double LumiPhotonCAL_XY = m_Params->get_double_param( "LumiPhotonCAL_XY" ) * cm;

  /*
  double LumiWin_X = m_Params->get_double_param( "LumiWin_X" ) * cm;
  double LumiWin_Y = m_Params->get_double_param( "LumiWin_Y" ) * cm;
  double LumiWin_Z = m_Params->get_double_param( "LumiWin_Z" ) * cm - enclosureCenter;
  double LumiWin_Tilt = m_Params->get_double_param( "LumiWin_Tilt" ) * rad;
  double LumiWin_Thickness = m_Params->get_double_param( "LumiWin_Thickness" ) * cm;
  double LumiWin_Height = m_Params->get_double_param( "LumiWin_Height" ) * cm;
  double LumiWin_Length = m_Params->get_double_param( "LumiWin_Length" ) * cm;
  std::string LumiWin_Material = m_Params->get_string_param( "LumiWin_Material" );
  
  double LumiMag_Z = m_Params->get_double_param( "LumiMag_Z" ) * cm - enclosureCenter;
  double LumiMag_inner = m_Params->get_double_param( "LumiMag_innerR" ) * cm;
  double LumiMag_outer = m_Params->get_double_param( "LumiMag_outerR" ) * cm;
  double LumiMag_DZ = m_Params->get_double_param( "LumiMag_DZ" ) * cm;
  double LumiMag_B = m_Params->get_double_param( "LumiMag_B" ) * tesla;
  std::string LumiMag_VesselMaterial = m_Params->get_string_param( "LumiMag_VesselMaterial" );

  //double LumiSpec_DZ = m_Params->get_double_param( "LumiSpec_DZ" ) * cm;
  //double LumiSpec_YS = m_Params->get_double_param( "LumiSpec_YS" ) * cm;
  //double LumiPhotonCAL_DZ = m_Params->get_double_param( "LumiPhotonCAL_DZ" ) * cm;

 
  //// Create G4 solid volumes: exit Window, B-field core + vessel

  G4Box *window = new G4Box("LumiWin", LumiWin_Length/2., LumiWin_Height/2., LumiWin_Thickness/2.);
  G4Tubs *core = new G4Tubs("LumiDipoleCore", 0., LumiMag_inner, LumiMag_DZ/2., 0., 360.*deg);
  G4Tubs *vessel = new G4Tubs("LumiDipoleVessel", LumiMag_inner, LumiMag_outer, LumiMag_DZ/2., 0., 360.*deg);
  
  // Create G4 logical volumes
  G4LogicalVolume *logical_window = new G4LogicalVolume( window, GetDetectorMaterial( LumiWin_Material ), "LumiWin");
  G4LogicalVolume *logical_core = new G4LogicalVolume( core, GetDetectorMaterial( "G4_Galactic" ), "LumiDipoleCore");
  G4LogicalVolume *logical_vessel = new G4LogicalVolume( vessel, GetDetectorMaterial( LumiMag_VesselMaterial ), "LumiDipoleVessel");
  
  G4VisAttributes *vis_window = new G4VisAttributes( G4Color(0, 1, 0, 0.5) );
  G4VisAttributes *vis_vessel = new G4VisAttributes( G4Color(0, 1, 0, 0.5) );
  vis_window->SetForceSolid( true );
  vis_vessel->SetForceSolid( true );
  
  logical_window->SetVisAttributes( vis_window );
  logical_core->SetVisAttributes( G4VisAttributes::GetInvisible() );
  logical_vessel->SetVisAttributes( vis_vessel );
  
  // Magnetic field for the core
  G4UniformMagField *field = new G4UniformMagField( G4ThreeVector( LumiMag_B, 0, 0 ) );
  G4FieldManager *fman = new G4FieldManager();
  fman->SetDetectorField( field );
  fman->CreateChordFinder( field );
  logical_core->SetFieldManager(fman, true);

  G4RotationMatrix *rot_win = new G4RotationMatrix( G4ThreeVector(0, 1, 0), LumiWin_Tilt ); //is typedef to CLHEP::HepRotation
  
  G4VPhysicalVolume *physical_window = new G4PVPlacement( rot_win, G4ThreeVector(LumiWin_X, LumiWin_Y, LumiWin_Z ), 
      logical_window, "LumiWin", logicWorld, 0, false, OverlapCheck());
  G4VPhysicalVolume *physical_core = new G4PVPlacement( 0, G4ThreeVector(0., 0., LumiMag_Z), 
      logical_core, "LumiDipoleCore", logicWorld, 0, false, OverlapCheck() );
  G4VPhysicalVolume *physical_vessel = new G4PVPlacement( 0, G4ThreeVector(0., 0., LumiMag_Z), 
      logical_vessel, "LumiDipoleVessel", logicWorld, 0, false, OverlapCheck() );

  // Add the physical volumes to appropriate container
  m_PassivePhysicalVolumesSet.insert( physical_window );
  m_PassivePhysicalVolumesSet.insert( physical_core );
  m_PassivePhysicalVolumesSet.insert( physical_vessel );

  // Add virtual layers for diagnostics
  AddVirtualLayer( "Virt_BeforeLumiDipole", G4TwoVector(LumiMag_outer, LumiMag_outer), G4ThreeVector(0., 0., LumiMag_Z + LumiMag_DZ/2. + 0.1*cm), logicWorld ); 
  AddVirtualLayer( "Virt_AfterLumiDipole", G4TwoVector(3*LumiMag_outer, 3*LumiMag_outer), G4ThreeVector(0., 0., LumiMag_Z - LumiMag_DZ/2. - 60*cm), logicWorld );

*/

   AddVirtualLayer( "Virt_UpperPhotonSpec", G4TwoVector(LumiSpec_XY, LumiSpec_XY), G4ThreeVector(0., LumiSpec_XY/2. + LumiPhotonCAL_XY/2. + 0.01*cm , LumiSpec_Z), logicWorld ); 
  AddVirtualLayer( "Virt_LowerPhotonSpec", G4TwoVector(LumiSpec_XY, LumiSpec_XY), G4ThreeVector(0., -(LumiSpec_XY/2. + LumiPhotonCAL_XY/2. + 0.01*cm) , LumiSpec_Z), logicWorld );
  AddVirtualLayer( "Virt_PhotonSpec", G4TwoVector(LumiPhotonCAL_XY, LumiPhotonCAL_XY), G4ThreeVector(0., 0., LumiPhotonCAL_Z), logicWorld ); 

    return;
}

//______________________________________________________________..
void EICG4LumiDetector::SetParametersFromFile()
{

	std::ifstream infile;
        std::string line;

        std::string paramFile = m_Params->get_string_param("parameter_file");   
	infile.open( paramFile );

	if( ! infile.is_open() ) 
	{
		std::cout << "ERROR in EICG4LumiDetector: Failed to open parameter file " << paramFile << std::endl;
		gSystem->Exit(1);
	}

	while( std::getline(infile, line) ) {

	    std::string name;
	    std::string value;

	    std::istringstream iss( line );

	    // skip comment lines
	    if( line.find("#") != std::string::npos ) { continue; }
            std::cout<<line<<endl;
	    if( !(iss >> name >> value) ) {
		std::cout << "Could not decode " << line << std::endl;
		gSystem->Exit(1);
	    }

            if( m_Params->exist_string_param( name ) ) {
                m_Params->set_string_param(name, value);
            }
            else if( m_Params->exist_double_param( name ) ) {
                m_Params->set_double_param(name, std::stod(value) );
            }
            else if( m_Params->exist_int_param( name ) ) {
                m_Params->set_int_param(name, std::stoi(value) );
            }
            else { 
                std::cout << "input parameter not recognized.  Exiting!" << std::endl;
                gSystem->Exit(1);
            }

	}
}



//_______________________________________________________________
void EICG4LumiDetector::Print(const std::string &what) const
{
  std::cout << "EICG4Lumi Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}

//______________________________________________________________

PHParameters *EICG4LumiDetector::getParams()
{
  return m_Params;
}

//_______________________________________________________________
void EICG4LumiDetector::AddVirtualLayer( std::string name, G4TwoVector size, G4ThreeVector pos, G4LogicalVolume *logicWorld ) 
{
  double virtPlaneDepth = 0.001 * cm;

  G4Box *virt = new G4Box(name, size.x()/2., size.y()/2., virtPlaneDepth);
  
  // Create G4 logical volumes
  G4LogicalVolume *logical_virt = new G4LogicalVolume( virt, GetDetectorMaterial( "G4_Galactic" ), name );

  G4VisAttributes *vis_virt = new G4VisAttributes( G4Color(1, 1, 1, 0.4) );
  vis_virt->SetForceSolid( true );
  logical_virt->SetVisAttributes( vis_virt );

  G4VPhysicalVolume *physical_virt = new G4PVPlacement( 0, pos, 
      logical_virt, name, logicWorld, 0, false, OverlapCheck());

  m_VirtualPhysicalVolumesMap.insert( {physical_virt, 1} );

}
