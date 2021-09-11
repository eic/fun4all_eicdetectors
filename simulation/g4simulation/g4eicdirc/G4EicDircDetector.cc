#include "G4EicDircDetector.h"

#include "G4EicDircDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Trd.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4Sphere.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Element.hh>
#include <Geant4/G4LogicalBorderSurface.hh>
#include <Geant4/G4LogicalSkinSurface.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RunManager.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4UImanager.hh>
#include <Geant4/G4VUserDetectorConstruction.hh>
#include <Geant4/G4ProcessManager.hh>

#include <cmath>
#include <iostream>  // for operator<<, endl, bas...

class G4VSolid;
class PHCompositeNode;

G4EicDircDetector::G4EicDircDetector(PHG4Subsystem *subsys,
                                     PHCompositeNode *Node,
                                     PHParameters *parameters,
                                     const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_DisplayAction(dynamic_cast<G4EicDircDisplayAction *>(subsys->GetDisplayAction()))
{
}

//_______________________________________________________________

int G4EicDircDetector::IsInDetector(G4VPhysicalVolume *volume) const
{ 
  std::map<G4VPhysicalVolume *, int>::const_iterator iter = m_PhysicalVolumes_active.find(volume);
  if(iter != m_PhysicalVolumes_active.end())
    {
      return iter->second;
    }

  return 0;
}


void G4EicDircDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  // ---------- DIRC supoort stucture ---------------

  //G4Material *Air = G4Material::GetMaterial("G4_AIR");

  // positions
  /*G4double rMin = m_Params->get_double_param("rMin");  // center location of Al support plate
  G4double det_height = 2.1 * cm;
  G4double place_z = m_Params->get_double_param("place_z");
  G4double detlength = m_Params->get_double_param("length");
  
  G4VSolid *dirc_envelope_solid = new G4Cons("dirc_envelope_solid",
                                             rMin - det_height / 2 - 2 * cm, rMin + det_height / 2 + 2 * cm,
                                             rMin - det_height / 2 - 2 * cm, rMin + det_height / 2 + 2 * cm,
                                             detlength / 2.0,
                                             0, 2 * M_PI);

  DetectorLog_Det = new G4LogicalVolume(dirc_envelope_solid, Air, name_base + "_Log");
 
  G4VPhysicalVolume* wDetectorLog_Det = new G4PVPlacement(0, G4ThreeVector(0, 0, place_z), DetectorLog_Det, name_base + "_Physical", logicWorld, false, 0, overlapcheck_sector); // FullEnvelope
  m_PhysicalVolumes_active[wDetectorLog_Det] = 20;

  // Single module with length based on readout (contains 14 LGADs [counting across both sides] in x-direction and 6 in z-direction)
  G4double baseplate_length = 43.1 * mm;
  //G4double baseplate_width = 56.5 * mm / 2;
  G4double segmentlength = 6 * baseplate_length;  //(detlength - 10 * cm) / 6;//m_Params->get_double_param("length");

  G4VSolid *sol_module_envelope = new G4Trd("sol_module_envelope",
                                            sin(M_PI / 12.) * rMin, sin(M_PI / 12.) * (rMin + det_height),
                                            segmentlength / 2, segmentlength / 2,
                                            det_height / 2);

  log_module_envelope = new G4LogicalVolume(sol_module_envelope, Air, "log_module_envelope");

  G4double cooling_plate_height = 6.35 * mm;
  G4double support_height = 7 * cm;

  // G4Material* mat_carbonfiber = new G4Material("CarbonFiberSupport", 1.44 * g / cm3, 1);
  // mat_carbonfiber->AddElement(G4Element::GetElement("C"), 1);
  // G4double density;  //z=mean number of protons;
  // G4int ncomponents;
  // carbon+epoxy material
  // G4Material *cfrp_intt = new G4Material("CFRP_INTT", density = 1.69 * g / cm3, ncomponents = 3);
  // cfrp_intt->AddElement(G4Element::GetElement("C"), 10);
  // cfrp_intt->AddElement(G4Element::GetElement("H"), 6);
  // cfrp_intt->AddElement(G4Element::GetElement("O"), 1);

  // SUPPORT STRUCTURES
  G4double support_width = 1 * mm;
  // build components of single segment here
  G4VSolid *Sol_End_Support = new G4Trd("Sol_End_Support",
                                        sin(M_PI / 12.) * (rMin - support_height * 0.9) - 2 * mm, sin(M_PI / 12.) * (rMin) -4 * mm,  // x1, x2
                                        support_width / 2, support_width / 2,                                                        // length
                                        support_height * 0.73 / 2);                                                                  // height

  Log_End_Support = new G4LogicalVolume(Sol_End_Support, G4Material::GetMaterial("G4_Fe"), "Log_End_Support_Raw");

  // place End side and back side support structure for the segment
  // RegisterPhysicalVolume( new G4PVPlacement(0, G4ThreeVector(0, segmentlength/2-support_width/2, -support_height/2), Log_End_Support,
  //                     "Front_Support_Physical", log_module_envelope, false, 0, overlapcheck_sector), false);
  // RegisterPhysicalVolume( new G4PVPlacement(0, G4ThreeVector(0, -segmentlength/2+support_width/2, -support_height/2), Log_End_Support,
  //                     "Back_Support_Physical", log_module_envelope, false, 0, overlapcheck_sector), false);


  // place longitudinal supports left, middle and right side of sector
  G4VSolid *Sol_Longitudinal_Support = new G4Trd("Sol_Longitudinal_Support",
                                                 support_width / 2, support_width / 2,                    // x1, x2
                                                 segmentlength / 2 - 1 * mm, segmentlength / 2 - 1 * mm,  // length
                                                 support_height * 0.73 / 2);                              // height

  Log_Longitudinal_Support = new G4LogicalVolume(Sol_Longitudinal_Support, G4Material::GetMaterial("G4_Fe"), "Log_Longitudinal_Support_Raw");

  G4RotationMatrix *supportrot = new G4RotationMatrix();
  supportrot->rotateY(-M_PI / 12.);
  if (rMin < 85 * cm)
    {
      G4VPhysicalVolume* wMother_Segment_Raw_Physical_Left = new G4PVPlacement(supportrot, G4ThreeVector(sin(M_PI / 12.) * (rMin - support_height / 2) - support_width / 2, 0, -support_height / 2), Log_Longitudinal_Support, "Mother_Segment_Raw_Physical_Left", log_module_envelope, false, 0, overlapcheck_sector); // Support
      m_PhysicalVolumes_active[wMother_Segment_Raw_Physical_Left] = 21;

      G4RotationMatrix *supportrot2 = new G4RotationMatrix();
      supportrot2->rotateY(M_PI / 12.);

      G4VPhysicalVolume* wMother_Segment_Raw_Physical_Right = new G4PVPlacement(supportrot2, G4ThreeVector(-sin(M_PI / 12.) * (rMin - support_height / 2) + support_width / 2, 0, -support_height / 2), Log_Longitudinal_Support, "Mother_Segment_Raw_Physical_Right", log_module_envelope, false, 0, overlapcheck_sector);
      m_PhysicalVolumes_active[wMother_Segment_Raw_Physical_Right] = 22;
    }

  G4double modulesep = 1 * mm;
  G4double moduleShift = -8 * mm;
  if (rMin < 85 * cm) moduleShift = -3 * mm;
  if (rMin < 66 * cm) moduleShift = -1 * mm;
  if (rMin < 55 * cm) moduleShift = 4 * mm;

  for (int isec = 0; isec < 12; isec++)
    {
      // if(isec!=3 && isec!=4)continue; // NOTE REMOVE
      // if(isec!=3)continue; // NOTE REMOVE
      G4RotationMatrix *motherrot = new G4RotationMatrix();
      motherrot->rotateX(M_PI / 2);
      motherrot->rotateY((isec - 3) * 2 * M_PI / 12.);
      // // central segments
      G4VPhysicalVolume* wlog_module_envelope = new G4PVPlacement(motherrot, G4ThreeVector((rMin - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), 0 * modulesep), log_module_envelope, "Mother_Segment_Raw_Physical_Center_" + std::to_string(isec), DetectorLog_Det, false, 0, overlapcheck_sector); // ModuleEnvelope
      
      m_PhysicalVolumes_active[wlog_module_envelope] = 23;
      
      for (int ilen = 1; ilen < ((detlength / 2 - segmentlength / 2) / segmentlength); ilen++)
	{
	  G4RotationMatrix *supfinalrot = new G4RotationMatrix();
	  // supfinalrot->rotateX(M_PI/2);
	  supfinalrot->rotateX(M_PI / 2);
	  supfinalrot->rotateY((isec - 3) * 2 * M_PI / 12.);
	  if (ilen == 2 || (ilen == 7))
	    {
	      if (rMin < 85 * cm)
		{		  
		  G4VPhysicalVolume* wFront_Support_Physical_1 = new G4PVPlacement(supfinalrot, G4ThreeVector((rMin - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * sin(isec* 2 * M_PI / 12.), ilen * segmentlength + segmentlength / 2 + ilen * modulesep), Log_End_Support, "Front_Support_Physical_1_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector);
		  m_PhysicalVolumes_active[wFront_Support_Physical_1] = 24;
		 
		  G4VPhysicalVolume* wFront_Support_Physical_2 = new G4PVPlacement(supfinalrot, G4ThreeVector((rMin - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * sin(isec* 2 * M_PI / 12.), -(ilen * segmentlength + segmentlength / 2 + ilen * modulesep)), Log_End_Support, "Front_Support_Physical_2_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector);

		  m_PhysicalVolumes_active[wFront_Support_Physical_2] = 25;
		    
		}    
	    }	
	    
	  // forward segments	  
	  G4VPhysicalVolume* wMother_Segment_Raw_Physical_Fwd = new G4PVPlacement(motherrot, G4ThreeVector((rMin - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), ilen * segmentlength + ilen * modulesep), log_module_envelope, "Mother_Segment_Raw_Physical_Fwd_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector);
	  
	  m_PhysicalVolumes_active[wMother_Segment_Raw_Physical_Fwd] = 26;


	  // backward segments	  
	  G4VPhysicalVolume* wMother_Segment_Raw_Physical_Bwd = new G4PVPlacement(motherrot, G4ThreeVector((rMin - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), -ilen * segmentlength - ilen * modulesep), log_module_envelope, "Mother_Segment_Raw_Physical_Bwd_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector);
	  
	  m_PhysicalVolumes_active[wMother_Segment_Raw_Physical_Bwd] = 27;

	}
    }

  // -------- INNER FRAME -----------------

  G4double rMin_inner = m_Params->get_double_param("rMin_inner");  // center location of Al support plate                                            
             
  G4VSolid *sol_module_envelope_inner = new G4Trd("sol_module_envelope_inner", sin(M_PI / 12.) * rMin_inner, sin(M_PI / 12.) * (rMin_inner + det_height), segmentlength / 2, segmentlength / 2, det_height / 2);

  log_module_envelope_inner = new G4LogicalVolume(sol_module_envelope_inner, Air, "log_module_envelope_inner");

  // SUPPORT STRUCTURES

  // build components of single segment here
  G4VSolid *Sol_End_Support_inner = new G4Trd("Sol_End_Support_inner", sin(M_PI / 12.) * (rMin_inner - support_height * 0.9) - 2 * mm, sin(M_PI / 12.) * (rMin_inner) -4 * mm,  support_width / 2, support_width / 2, support_height * 0.73 / 2);                                                           

  Log_End_Support_inner = new G4LogicalVolume(Sol_End_Support_inner, G4Material::GetMaterial("G4_Fe"), "Log_End_Support_Raw_inner");
  
  if (rMin_inner < 85 * cm)
    {
      G4VPhysicalVolume* wMother_Segment_Raw_Physical_Left_inner = new G4PVPlacement(supportrot, G4ThreeVector(sin(M_PI / 12.) * (rMin_inner - support_height / 2) - support_width / 2, 0, -support_height / 2), Log_Longitudinal_Support, "Mother_Segment_Raw_Physical_Left_inner", log_module_envelope_inner, false, 0, overlapcheck_sector); // Support
      m_PhysicalVolumes_active[wMother_Segment_Raw_Physical_Left_inner] = 28;

      G4RotationMatrix *supportrot2 = new G4RotationMatrix();
      supportrot2->rotateY(M_PI / 12.);

      G4VPhysicalVolume* wMother_Segment_Raw_Physical_Right_inner = new G4PVPlacement(supportrot2, G4ThreeVector(-sin(M_PI / 12.) * (rMin_inner - support_height / 2) + support_width / 2, 0, -support_height / 2), Log_Longitudinal_Support, "Mother_Segment_Raw_Physical_Right_inner", log_module_envelope_inner, false, 0, overlapcheck_sector);
      m_PhysicalVolumes_active[wMother_Segment_Raw_Physical_Right_inner] = 29;
    }

  //G4double modulesep = 1 * mm;
  //G4double moduleShift = -8 * mm;
  if (rMin_inner < 85 * cm) moduleShift = -3 * mm;
  //if (rMin < 66 * cm) moduleShift = -1 * mm;
  //if (rMin < 55 * cm) moduleShift = 4 * mm;

  for (int isec = 0; isec < 12; isec++)
    {
      // if(isec!=3 && isec!=4)continue; // NOTE REMOVE
      // if(isec!=3)continue; // NOTE REMOVE
      G4RotationMatrix *motherrot = new G4RotationMatrix();
      motherrot->rotateX(M_PI / 2);
      motherrot->rotateY((isec - 3) * 2 * M_PI / 12.);
      // // central segments
      G4VPhysicalVolume* wlog_module_envelope_inner = new G4PVPlacement(motherrot, G4ThreeVector((rMin_inner - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin_inner - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), 0 * modulesep), log_module_envelope_inner, "Mother_Segment_Raw_Physical_Center_inner_" + std::to_string(isec), DetectorLog_Det, false, 0, overlapcheck_sector); // ModuleEnvelope
      
      m_PhysicalVolumes_active[wlog_module_envelope_inner] = 30;
      
      for (int ilen = 1; ilen < ((detlength / 2 - segmentlength / 2) / segmentlength); ilen++)
	{
	  G4RotationMatrix *supfinalrot = new G4RotationMatrix();
	  // supfinalrot->rotateX(M_PI/2);
	  supfinalrot->rotateX(M_PI / 2);
	  supfinalrot->rotateY((isec - 3) * 2 * M_PI / 12.);
	  if (ilen == 2 || (ilen == 7))
	    {
	      if (rMin_inner < 85 * cm)
		{  
		  G4VPhysicalVolume* wFront_Support_Physical_inner_1 = new G4PVPlacement(supfinalrot, G4ThreeVector((rMin_inner - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin_inner - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * sin(isec* 2 * M_PI / 12.), ilen * segmentlength + segmentlength / 2 + ilen * modulesep), Log_End_Support_inner, "Front_Support_Physical_inner_1_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector);
		  m_PhysicalVolumes_active[wFront_Support_Physical_inner_1] = 31;
		   
		  G4VPhysicalVolume* wFront_Support_Physical_inner_2 = new G4PVPlacement(supfinalrot, G4ThreeVector((rMin_inner - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin_inner - support_height / 2 - det_height / 2 - cooling_plate_height / 2 + moduleShift) * sin(isec* 2 * M_PI / 12.), -(ilen * segmentlength + segmentlength / 2 + ilen * modulesep)), Log_End_Support_inner, "Front_Support_Physical_inner_2_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector);

		  m_PhysicalVolumes_active[wFront_Support_Physical_inner_2] = 32;
		      
		}    
	    }
	      
	  // forward segments  
	  G4VPhysicalVolume* wMother_Segment_Raw_Physical_Fwd_inner = new G4PVPlacement(motherrot, G4ThreeVector((rMin_inner - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin_inner - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), ilen * segmentlength + ilen * modulesep), log_module_envelope_inner, "Mother_Segment_Raw_Physical_Fwd_inner_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector);
	    
	  m_PhysicalVolumes_active[wMother_Segment_Raw_Physical_Fwd_inner] = 33;


	  // backward segments  
	  G4VPhysicalVolume* wMother_Segment_Raw_Physical_Bwd_inner = new G4PVPlacement(motherrot, G4ThreeVector((rMin_inner - det_height / 2 + moduleShift) * cos(isec * 2 * M_PI / 12.), (rMin_inner - det_height / 2 + moduleShift) * sin(isec * 2 * M_PI / 12.), -ilen * segmentlength - ilen * modulesep), log_module_envelope_inner, "Mother_Segment_Raw_Physical_Bwd_inner_" + std::to_string(isec) + "_" + std::to_string(ilen), DetectorLog_Det, false, 0, overlapcheck_sector);
	    
	  m_PhysicalVolumes_active[wMother_Segment_Raw_Physical_Bwd_inner] = 34;

	}
	}*/


  // -------------- DIRC ----------------------

  fGeomType = m_Params->get_int_param("Geom_type");
  fLensId = m_Params->get_int_param("Lens_id"); 
  fNBar = m_Params->get_double_param("NBars");

  std::cout << "Nbars = " << fNBar << std::endl;

  fNRow = m_Params->get_int_param("MCP_rows");
  fNCol = m_Params->get_int_param("MCP_columns");

  fPrizm[0] = m_Params->get_double_param("Prizm_width");
  fPrizm[1] = m_Params->get_double_param("Prizm_length"); 
  fPrizm[3] = m_Params->get_double_param("Prizm_height_at_lens");
  fPrizm[2]= fPrizm[3] + (fPrizm[1]*tan(32*deg));
  
  fBar[0] = fBarL[0] = fBarS[0] = m_Params->get_double_param("Bar_thickness"); 
  fBar[1] = fBarL[1] = fBarS[1] = m_Params->get_double_param("Bar_width"); 
  fBarL[2] = m_Params->get_double_param("BarL_length");
  fBarS[2] = m_Params->get_double_param("BarS_length");
  
  fNBoxes = m_Params->get_int_param("NBoxes");
  fRadius = m_Params->get_double_param("Radius");
  zshift = m_Params->get_double_param("z_shift");

  //-------------------

  fBarsGap = 0.15;
  
  std::cout << "Nbarboxes = "<< fNBoxes << std::endl;
  std::cout << "barrel radius = " << fRadius << " mm" << std::endl;
  
  fMirror[0] = m_Params->get_double_param("Mirror_height"); 
  fMirror[1] = fPrizm[0]; fMirror[2] =1;
  
  fMcpTotal[0] = fMcpTotal[1] = 53+4; fMcpTotal[2]=1;
  fMcpActive[0] = fMcpActive[1] = 53; fMcpActive[2]=1;
  fLens[0] = fLens[1] = 40; fLens[2]=10;
  
  fBoxWidth = fPrizm[0];

  if(fGeomType == 1)  fNBoxes = 1;

  fFd[0] = fBoxWidth; fFd[1]=fPrizm[2]; fFd[2] =1;

  
  if(fLensId == 0 || fLensId==10){
    fLens[2] =0;
  }
  if(fLensId == 2){
    fLens[0] = fPrizm[3]; fLens[1] = 175; fLens[2]=14.4;
  }
  if(fLensId == 3){
    fLens[0] = fPrizm[3]; fLens[1] = fPrizm[0]/fNBar; fLens[2]=15;
    if(fNBar==1)  fLens[1] = fPrizm[0]/11.;
  }

  if(fLensId == 6){
    fLens[0] = fPrizm[3]; fLens[1] = fPrizm[0]; fLens[2]=12;
  }  

  fPrtRot = new G4RotationMatrix();
  
  //-------------------------
  DefineMaterials();

  // ------------- Volumes --------------

  double gluethickness = 0.05;
  int nparts = m_Params->get_int_param("Bar_pieces");

  double dirclength = fBarL[2]*(nparts-1) + fBarS[2] + gluethickness*nparts;

  // The DIRC
  G4Box* gDirc = new G4Box("gDirc",210,195,0.5*dirclength+350);
  lDirc = new G4LogicalVolume(gDirc,defaultMaterial,"lDirc",0,0,0);

  G4Box* gFd = new G4Box("gFd",0.5*fFd[1],0.5*fFd[0],0.5*fFd[2]);
  lFd = new G4LogicalVolume(gFd,defaultMaterial,"lFd",0,0,0);

  double dphi = 360*deg/(double)fNBoxes;  
  G4VPhysicalVolume* phy;


  // LUT
  /*phy = new G4PVPlacement(0,G4ThreeVector(0,0,0),lDirc,"wDirc",logicWorld,false,0);
    m_PhysicalVolumesSet.insert(phy);
  */ 
    for(int i=0; i<fNBoxes; i++){
      double tphi = dphi*i; 
      double dx = fRadius * cos(tphi);
      double dy = fRadius * sin(tphi);

      G4RotationMatrix *tRot = new G4RotationMatrix();
      tRot->rotateZ(-tphi);     
      phy = new G4PVPlacement(tRot,G4ThreeVector(dx,dy,zshift),lDirc,"wDirc",logicWorld,false,i,OverlapCheck());
      m_PhysicalVolumes_active[phy] = 1;
    }
  
  // The Bar
  G4Box *gBarL = new G4Box("gBarL", fBarL[0]/2., fBarL[1]/2., fBarL[2]/2.);
  lBarL = new G4LogicalVolume(gBarL,BarMaterial,"lBarL",0,0,0);

  G4Box *gBarS = new G4Box("gBarS", fBarS[0]/2., fBarS[1]/2., fBarS[2]/2.);
  lBarS = new G4LogicalVolume(gBarS,BarMaterial,"lBarS",0,0,0);
    
  // Glue
  G4Box *gGlue = new G4Box("gGlue",fBar[0]/2.,fBar[1]/2.,0.5*gluethickness);
  lGlue = new G4LogicalVolume(gGlue,epotekMaterial,"lGlue",0,0,0);

  int id=0; 

  for(int i=0; i<fNBar; i++){
    double shifty = i*(fBar[1]+fBarsGap)- 0.5*fBoxWidth + fBar[1]/2.; 
    for(int j=0; j<nparts; j++){
      if(j < (nparts-1))
	{
	  double z = 0.5*dirclength - 0.5*fBarL[2] - (fBarL[2]+gluethickness)*j;
	  pDirc[i] = new G4PVPlacement(0,G4ThreeVector(0,shifty,z),lBarL,"wBar",lDirc,false,id,OverlapCheck());
	  m_PhysicalVolumes_active[pDirc[i]] = 3;
	  wGlue = new G4PVPlacement(0,G4ThreeVector(0,shifty,z-0.5*(fBarL[2]+gluethickness)),lGlue,"wGlue", lDirc,false,id,OverlapCheck());
	  m_PhysicalVolumes_active[wGlue] = 4;
	}
      else
	{
	  double z = 0.5*dirclength - (nparts-1)*fBarL[2] - 0.5*fBarS[2] - gluethickness*j;
	  pDirc[i] = new G4PVPlacement(0,G4ThreeVector(0,shifty,z),lBarS,"wBar",lDirc,false,id,OverlapCheck());
	  m_PhysicalVolumes_active[pDirc[i]] = 3;
	  wGlue = new G4PVPlacement(0,G4ThreeVector(0,shifty,z-0.5*(fBarS[2]+gluethickness)),lGlue,"wGlue", lDirc,false,id,OverlapCheck());
	  m_PhysicalVolumes_active[wGlue] = 4;
	}
      id++;
    }
  }

  
  // The Mirror
  G4Box* gMirror = new G4Box("gMirror",fMirror[0]/2.,fMirror[1]/2.,fMirror[2]/2.);
  lMirror = new G4LogicalVolume(gMirror,MirrorMaterial,"lMirror",0,0,0);
  wMirror =new G4PVPlacement(0,G4ThreeVector(0,0,0.5*dirclength+fMirror[2]/2.),lMirror,"wMirror", lDirc,false,0,OverlapCheck());
  m_PhysicalVolumes_active[wMirror] = 5;  

  // The Lenses
  if(fLensId == 2){ // 2-layer lens
    double lensrad = 70;
    double lensMinThikness = 2;
    G4Box* gfbox = new G4Box("Fbox",fLens[0]/2.,fLens[1]/2.,fLens[2]/2.);
    G4ThreeVector zTrans(0, 0, -lensrad+fLens[2]/2.-lensMinThikness);

    G4Sphere* gsphere = new G4Sphere("Sphere",0,70,0.*deg,360.*deg,0.*deg,380.*deg);
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Fbox*Sphere", gfbox, gsphere,new G4RotationMatrix(),zTrans); 
    G4SubtractionSolid* gLens2 = new G4SubtractionSolid("Fbox-Sphere", gfbox, gsphere,new G4RotationMatrix(),zTrans);

    lLens1 = new G4LogicalVolume(gLens1,Nlak33aMaterial,"lLens1",0,0,0); //Nlak33aMaterial  
    lLens2 = new G4LogicalVolume(gLens2,BarMaterial,"lLens2",0,0,0);
  }

  if(fLensId == 3){ // 3-component spherical lens
    double lensMinThikness = 2; 
  
    double r1 = 47.8;//0; 
    double r2 = 29.1;//0; 
  
    //r1 = (r1==0)? 47.8: r1;
    //r2 = (r2==0)? 29.1: r2;
    // r1=80;
    // r2=35;

    double cr2 = sqrt(fLens[1]*fLens[1]/4.+fBar[0]*fBar[0]/4.);
    if(cr2 > r2) {
      std::cout<<"bad lens"<<std::endl;
      cr2 = r2;
    }
    fLens[2] = (2*lensMinThikness+r2-sqrt(r2*r2-cr2*cr2)+lensMinThikness);
    std::cout << "lens thickness ="<< fLens[2] << " mm" << std::endl;
    
    G4ThreeVector zTrans1(0, 0, -r1-fLens[2]/2.+r1-sqrt(r1*r1-cr2*cr2) +lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -r2-fLens[2]/2.+r2-sqrt(r2*r2-cr2*cr2) +lensMinThikness*2);

    G4Box* gfbox = new G4Box("Fbox",fLens[0]/2.,fLens[1]/2.,fLens[2]/2.);
    G4Tubs* gfstube = new G4Tubs("ftube",0,cr2,fLens[2]/2.,0,360*deg);

    G4Sphere* gsphere1 = new G4Sphere("Sphere1",0,r1,0,360*deg,0,360*deg);
    G4Sphere* gsphere2 = new G4Sphere("Sphere2",0,r2,0,360*deg,0,360*deg);

    G4IntersectionSolid* gbbox = new G4IntersectionSolid("bbox", gfbox, gfbox,new G4RotationMatrix(),G4ThreeVector(0,0,-lensMinThikness*2)); 
    G4IntersectionSolid* gsbox = new G4IntersectionSolid("sbox", gfstube, gfbox,new G4RotationMatrix(),G4ThreeVector(0,0,lensMinThikness*2)); 


    G4UnionSolid* gubox = new G4UnionSolid("unionbox", gbbox, gsbox,new G4RotationMatrix(),G4ThreeVector(0,0,0)); 

    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Lens1", gubox, gsphere1,new G4RotationMatrix(),-zTrans1); 
    G4SubtractionSolid*  gLenst = new G4SubtractionSolid("temp", gubox, gsphere1, new G4RotationMatrix(),-zTrans1);

    G4IntersectionSolid* gLens2 = new G4IntersectionSolid("Lens2", gLenst, gsphere2, new G4RotationMatrix(),-zTrans2);
    G4SubtractionSolid*  gLens3 = new G4SubtractionSolid("Lens3", gLenst, gsphere2,new G4RotationMatrix(),-zTrans2);
    
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0);
    lLens2 = new G4LogicalVolume(gLens2,Nlak33aMaterial,"lLens2",0,0,0); //Nlak33aMaterial //PbF2Material //SapphireMaterial
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }

  if(fLensId == 6){ // 3-component cylindrical lens
    double lensMinThikness = 2.0;

    double r1 = 33;//0; 
    double r2 = 24;//0; 

    lensMinThikness = 2;
    double layer12 = lensMinThikness*2;

    // r1 = (r1==0)? 27.45: r1;
    // r2 = (r2==0)? 20.02: r2;

    //r1 = (r1==0)? 33: r1;
    //r2 = (r2==0)? 24: r2;
    double shight = 25;

    G4ThreeVector zTrans1(0, 0, -r1-fLens[2]/2.+r1-sqrt(r1*r1-shight/2.*shight/2.) +lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -r2-fLens[2]/2.+r2-sqrt(r2*r2-shight/2.*shight/2.) +layer12);

    G4Box* gfbox = new G4Box("fbox",0.5*fLens[0],0.5*fLens[1],0.5*fLens[2]);
    G4Box* gcbox = new G4Box("cbox",0.5*fLens[0],0.5*fLens[1]+1,0.5*fLens[2]);
    G4ThreeVector tTrans1( 0.5*(fLens[0]+shight),0,-fLens[2]+layer12);
    G4ThreeVector tTrans0(-0.5*(fLens[0]+shight),0,-fLens[2]+layer12);
    G4SubtractionSolid*  tubox = new G4SubtractionSolid("tubox", gfbox, gcbox,new G4RotationMatrix(),tTrans1);
    G4SubtractionSolid*  gubox = new G4SubtractionSolid("gubox", tubox, gcbox,new G4RotationMatrix(),tTrans0);

    G4Tubs* gcylinder1  = new G4Tubs("Cylinder1",0,r1,0.5*fLens[1],0*deg,360*deg);
    G4Tubs* gcylinder2  = new G4Tubs("Cylinder2",0,r2,0.5*fLens[1]-0.5,0*deg,360*deg);
    G4Tubs* gcylinder1c = new G4Tubs("Cylinder1c",0,r1,0.5*fLens[1]+0.5,0*deg,360*deg);
    G4Tubs* gcylinder2c = new G4Tubs("Cylinder2c",0,r2,0.5*fLens[1]+0.5,0*deg,360*deg);
    G4RotationMatrix* xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI/2.*rad);

    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Lens1", gubox, gcylinder1,xRot,zTrans1);
    G4SubtractionSolid*  gLenst = new G4SubtractionSolid("temp", gubox, gcylinder1c, xRot,zTrans1);

    G4IntersectionSolid* gLens2 = new G4IntersectionSolid("Lens2", gLenst, gcylinder2, xRot,zTrans2);
    G4SubtractionSolid*  gLens3 = new G4SubtractionSolid("Lens3", gLenst, gcylinder2c, xRot,zTrans2);

    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0);
    lLens2 = new G4LogicalVolume(gLens2,Nlak33aMaterial,"lLens2",0,0,0);
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }  

  if(fLensId ==100){
    fLens[2]=200;
    G4Box *gLens3 = new G4Box("gLens1",fBar[0]/2.,0.5*fBoxWidth,100);
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }
      
  if(fLensId != 0 && fLensId != 10){
    double shifth = -0.5*(dirclength+fLens[2]);
    if(fLensId==100){
      phy = new G4PVPlacement(0,G4ThreeVector(0,0,shifth),lLens3,"wLens3", lDirc,false,0,OverlapCheck());

    }else if(fNBar==1 && fLensId==3){
      for(int i=0; i<11; i++){
	double shifty = i*fLens[1] - 0.5*(fBoxWidth - fLens[1]);
	phy = new G4PVPlacement(0,G4ThreeVector(0,shifty,shifth),lLens1,"wLens1", lDirc,false,i,OverlapCheck());
	phy = new G4PVPlacement(0,G4ThreeVector(0,shifty,shifth),lLens2,"wLens2", lDirc,false,i,OverlapCheck());
	phy = new G4PVPlacement(0,G4ThreeVector(0,shifty,shifth),lLens3,"wLens3", lDirc,false,i,OverlapCheck());

      }
    }else{
      for(int i=0; i<fNBar; i++){
	double shifty = i*fLens[1] - 0.5*(fBoxWidth - fLens[1]);
	if(fLensId!=6){
	  phy = new G4PVPlacement(0,G4ThreeVector(0,shifty,shifth),lLens1,"wLens1", lDirc,false,i,OverlapCheck());
	  m_PhysicalVolumes_active[phy] = 6;

	  phy = new G4PVPlacement(0,G4ThreeVector(0,shifty,shifth),lLens2,"wLens2", lDirc,false,i,OverlapCheck());
	  m_PhysicalVolumes_active[phy] = 7;

	  if(fLensId == 3)  phy = new G4PVPlacement(0,G4ThreeVector(0,shifty,shifth),lLens3,"wLens3", lDirc,false,i,OverlapCheck());
	  m_PhysicalVolumes_active[phy] = 8; 
	}
      }
      if(fLensId==6){
	double sh=0;
	phy = new G4PVPlacement(0,G4ThreeVector(0,0,shifth-sh),lLens1,"wLens1", lDirc,false,0,OverlapCheck());
	phy = new G4PVPlacement(0,G4ThreeVector(0,0,shifth-sh),lLens2,"wLens2", lDirc,false,0,OverlapCheck());
	phy = new G4PVPlacement(0,G4ThreeVector(0,0,shifth-sh),lLens3,"wLens3", lDirc,false,0,OverlapCheck());

      }
    }
  }

  // The Prizm
  G4Trap* gPrizm = new G4Trap("gPrizm",fPrizm[0],fPrizm[1],fPrizm[2],fPrizm[3]);
  lPrizm = new G4LogicalVolume(gPrizm, BarMaterial,"lPrizm",0,0,0);

  G4RotationMatrix* xRot = new G4RotationMatrix();
  xRot->rotateX(-M_PI/2.*rad);


  G4RotationMatrix *fdrot = new G4RotationMatrix();
  double evshiftz = 0.5*dirclength+fPrizm[1]+fMcpActive[2]/2.+fLens[2];
  double evshiftx = 0;

  fPrismShift = G4ThreeVector((fPrizm[2]+fPrizm[3])/4.-0.5*fPrizm[3],0,-(0.5*(dirclength+fPrizm[1])+fLens[2]));
  phy = new G4PVPlacement(xRot,fPrismShift,lPrizm,"wPrizm", lDirc,false,0,OverlapCheck());
  m_PhysicalVolumes_active[phy] = 9;

  phy = new G4PVPlacement(fdrot,G4ThreeVector(0.5*fFd[1]-0.5*fPrizm[3]-evshiftx,0,-evshiftz),lFd,"wFd", lDirc,false,0,OverlapCheck());
  m_PhysicalVolumes_active[phy] = 2;
  
  // MCP --

  G4Box* gMcp = new G4Box("gMcp",fMcpTotal[0]/2.,fMcpTotal[1]/2.,fMcpTotal[2]/2.);
  lMcp = new G4LogicalVolume(gMcp,BarMaterial,"lMcp",0,0,0);

  fNpix1 = 16;
  fNpix2 = 16;

  std::cout<<"fNpix1="<<fNpix1 << " fNpix2="<<fNpix2 <<std::endl;
        
  // The MCP Pixel
  G4Box* gPixel = new G4Box("gPixel",0.5*fMcpActive[0]/fNpix1,0.5*fMcpActive[1]/fNpix2,fMcpActive[2]/16.);
  lPixel = new G4LogicalVolume(gPixel,BarMaterial,"lPixel",0,0,0);

  for(int i=0; i<fNpix2; i++){
    for(int j=0; j<fNpix1; j++){
      double shiftx = i*(fMcpActive[0]/fNpix1)-fMcpActive[0]/2.+0.5*fMcpActive[0]/fNpix1;
      double shifty = j*(fMcpActive[1]/fNpix2)-fMcpActive[1]/2.+0.5*fMcpActive[1]/fNpix2;
      phy = new G4PVPlacement(0,G4ThreeVector(shiftx,shifty,0),lPixel,"wPixel", lMcp,false,fNpix2*i+j,OverlapCheck());      
      m_PhysicalVolumes_active[phy] = 11;
    }
  }
 
  int mcpId = 0;
  double gapx = (fPrizm[2]-fNCol*fMcpTotal[0])/(double)(fNCol+1);
  double gapy = (fBoxWidth-fNRow*fMcpTotal[1])/(double)(fNRow+1);
  for(int i=0; i<fNCol; i++){
    for(int j=0; j<fNRow; j++){

      double shiftx = i*(fMcpTotal[0]+gapx)-0.5*(fFd[1]-fMcpTotal[0])+gapx;
      double shifty = j*(fMcpTotal[1]+gapy)-0.5*(fBoxWidth-fMcpTotal[1])+gapy;
      //phy = new G4PVPlacement(0,G4ThreeVector(shiftx,shifty,0),lMcp,"wMcp", lFd,false,mcpId++);
      phy = new G4PVPlacement(0,G4ThreeVector(shiftx,shifty,0),lMcp,"wMcp", lFd,false,mcpId,OverlapCheck()); 
      m_PhysicalVolumes_active[phy] = 10;
      mcpId++;
    }
  }


  {
    const int num = 36; 
    double WaveLength[num];
    double PhotonEnergy[num];
    double PMTReflectivity[num];
    double EfficiencyMirrors[num];
    const double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
    for(int i=0;i<num;i++){
      WaveLength[i]= (300 +i*10)*nanometer;
      PhotonEnergy[num-(i+1)]= LambdaE/WaveLength[i];
      PMTReflectivity[i]=0.;
      EfficiencyMirrors[i]=0; 
    }

    /***************** QUANTUM EFFICIENCY OF BURLE AND HAMAMTSU PMT'S *****/

    //ideal pmt quantum efficiency
    double QuantumEfficiencyIdial[num]=
      {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
       1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
       1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
       1.0,1.0,1.0,1.0,1.0,1.0};

    // Burle PMT's 
    double QuantumEfficiencyB[num] =
      {0.,0.001,0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06,
       0.07,0.09,0.1,0.13,0.15,0.17,0.2,0.24,0.26,0.28,0.282,0.284,0.286,
       0.288,0.29,0.28,0.26,0.24,0.22,0.20,0.18,0.15,0.13,0.12,0.10};
  
    //hamamatsu pmt quantum efficiency
    double QuantumEfficiencyPMT[num]=
      {0.001,0.002,0.004,0.007,0.011,0.015,0.020,0.026,0.033,0.040,0.045,
       0.056,0.067,0.085,0.109,0.129,0.138,0.147,0.158,0.170,
       0.181,0.188,0.196,0.203,0.206,0.212,0.218,0.219,0.225,0.230,
       0.228,0.222,0.217,0.210,0.199,0.177};

    // these quantum efficiencies have to be multiplied by geometry
    //   efficiency of given PMT's
    //   for Hamamatsu by factor 0.7
    //   for Burle by factor 0.45 
    for(int k=0;k<36;k++){
      QuantumEfficiencyB[k] =  QuantumEfficiencyB[k] * 0.45 ;
      QuantumEfficiencyPMT[k] =  QuantumEfficiencyPMT[k] *.7;
    }
 
    
    /* define quantum efficiency for burle PMT's - the same efficiency is 
       assigned to pads and also to slots !!!! */
    
    //burle pmt - bigger slots => logicPad
    G4MaterialPropertiesTable* PhotocatodBurleMPT = new G4MaterialPropertiesTable();
    PhotocatodBurleMPT->AddProperty("EFFICIENCY",  PhotonEnergy,QuantumEfficiencyB,num);
    PhotocatodBurleMPT->AddProperty("REFLECTIVITY",PhotonEnergy,PMTReflectivity,num);

 
    G4OpticalSurface* BurlePMTOpSurface= 
      new G4OpticalSurface("BurlePMTOpSurface",glisur,polished,
			   dielectric_metal);
    BurlePMTOpSurface->SetMaterialPropertiesTable(PhotocatodBurleMPT);

    
    /* hamamatsu pmt's - smaller slots => quantum efficiency again
       assign to slot and pad */
  
    fQuantumEfficiency = QuantumEfficiencyIdial;
    G4MaterialPropertiesTable* PhotocatodHamamatsuMPT = new G4MaterialPropertiesTable();
    PhotocatodHamamatsuMPT->AddProperty("EFFICIENCY",  PhotonEnergy,fQuantumEfficiency,num);
    PhotocatodHamamatsuMPT->AddProperty("REFLECTIVITY",PhotonEnergy,PMTReflectivity,num);

    G4OpticalSurface* HamamatsuPMTOpSurface= 
      new G4OpticalSurface("HamamatsuPMTOpSurface",glisur,polished,dielectric_metal);
    HamamatsuPMTOpSurface->SetMaterialPropertiesTable(PhotocatodHamamatsuMPT);

    // // assignment to pad
    // if(hamamatsu8500)
    new G4LogicalSkinSurface("HamamatsuPMTSurface",lPixel,HamamatsuPMTOpSurface);

    // Mirror
    G4OpticalSurface* MirrorOpSurface= 
      new G4OpticalSurface("MirrorOpSurface",glisur,polished,dielectric_metal);
  
    double ReflectivityMirrorBar[num]={
      0.8755,0.882,0.889,0.895,0.9,0.9025,0.91,0.913,0.9165,0.92,0.923,
      0.9245,0.9285,0.932,0.934,0.935,0.936,0.9385,0.9395,0.94,
      0.9405,0.9405,0.9405,0.9405,0.94,0.9385,0.936,0.934,
      0.931,0.9295,0.928,0.928,0.921,0.92,0.927,0.9215};

    G4MaterialPropertiesTable *MirrorMPT = new G4MaterialPropertiesTable();
    MirrorMPT->AddProperty("REFLECTIVITY", PhotonEnergy, ReflectivityMirrorBar, num);
    MirrorMPT->AddProperty("EFFICIENCY", PhotonEnergy, EfficiencyMirrors,   num);
  
    MirrorOpSurface->SetMaterialPropertiesTable(MirrorMPT);
    new G4LogicalSkinSurface("MirrorSurface", lMirror,MirrorOpSurface);
  }
  
  SetVisualization();

  /*PrtOpBoundaryProcess *fBoundaryProcess = new PrtOpBoundaryProcess();
  G4ProcessManager *pmanager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  pmanager->AddDiscreteProcess(fBoundaryProcess);
  pmanager->SetProcessOrderingToFirst(fBoundaryProcess, idxPostStep);
  */

}

void G4EicDircDetector::DefineMaterials()
{
  G4String symbol;             //a=mass of a mole;
  double a, z, density;      //z=mean number of protons;  

  int ncomponents, natoms;
  double fractionmass;

  // define Elements
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silicon" ,symbol="Si", z= 14., a= 28.09*g/mole);
  G4Element* Al = new G4Element("Aluminum",symbol="Al", z=13., a=26.98*g/mole);


  // quartz material = SiO2
  G4Material* SiO2 = new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);

  Nlak33aMaterial  = new G4Material("Nlak33a",density= 4.220*g/cm3, ncomponents=2);
  Nlak33aMaterial->AddElement(Si, natoms=1);
  Nlak33aMaterial->AddElement(O , natoms=2);

  PbF2Material  = new G4Material("PbF2",density= 3.97*g/cm3, ncomponents=2);
  PbF2Material->AddElement(Si, natoms=1);
  PbF2Material->AddElement(O , natoms=2);

  SapphireMaterial  = new G4Material("Sapphire",density= 4.220*g/cm3, ncomponents=2);
  SapphireMaterial->AddElement(Al, natoms=2);
  SapphireMaterial->AddElement(O , natoms=3);
  
  G4Material* Vacuum = new G4Material("interGalactic", 1., 1.008*g/mole, 
				      1.e-25*g/cm3, kStateGas, 
				      2.73*kelvin, 3.e-18*pascal);
  G4Material* Air = new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Aluminum = new G4Material("Aluminum",density=2.7*g/cm3,ncomponents=1);
  Aluminum->AddElement(Al,fractionmass=1.0);

  G4Material* KamLandOil = new G4Material("KamLandOil",density=0.914*g/cm3,ncomponents=2);
  KamLandOil->AddElement(C,natoms=12);
  KamLandOil->AddElement(H,natoms=26);

  G4Material* CarbonFiber = new G4Material("CarbonFiber", density=0.145*g/cm3, ncomponents=1);
  CarbonFiber->AddElement(C,natoms=1);
  

  /* as I don't know the exact material composition,
     I will use Epoxyd material composition and add
     the optical property of Epotek to this material */

  G4Material* Epotek = new G4Material("Epotek",density=1.2*g/cm3,ncomponents=3);

  Epotek->AddElement(C,natoms=3);
  Epotek->AddElement(H,natoms=5);
  Epotek->AddElement(O,natoms=2);


  // assign main materials
  if(fGeomType < 10) defaultMaterial = Vacuum;
  else defaultMaterial = Air; //Vacuum // material of world
  frontMaterial =  CarbonFiber; 
  BarMaterial = SiO2; // material of all Bars, Quartz and Window
  OilMaterial = KamLandOil; // material of volume 1,2,3,4
  MirrorMaterial = Aluminum; // mirror material
  epotekMaterial = Epotek; // Epotek material - glue between bars


  // ------------ Generate & Add Material Properties Table ------------

  static const double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
  const int num = 36;
  double WaveLength[num];
  //double Absorption[num]; // default value for absorption
  double AirAbsorption[num]; // absorption value for air
  double AirRefractiveIndex[num]; // air refractive index
  double PhotonEnergy[num]; // energy of photons which correspond to the given 
  // refractive or absoprtion values

  double PhotonEnergyNlak33a[76] = {1,1.2511,1.26386,1.27687,1.29016,1.30372,1.31758,1.33173,1.34619,1.36097,1.37607,1.39152,1.40731,1.42347,1.44,1.45692,1.47425,1.49199,1.51016,1.52878,1.54787,1.56744,1.58751,1.6081,1.62923,1.65092,1.6732,1.69609,1.71961,1.7438,1.76868,1.79427,1.82062,1.84775,1.87571,1.90452,1.93423,1.96488,1.99652,2.0292,2.06296,2.09787,2.13398,2.17135,2.21006,2.25017,2.29176,2.33492,2.37973,2.42631,2.47473,2.52514,2.57763,2.63236,2.68946,2.7491,2.81143,2.87666,2.94499,3.01665,3.09187,3.17095,3.25418,3.34189,3.43446,3.53231,3.6359,3.74575,3.86244,3.98663,4.11908,4.26062,4.41225,4.57506,4.75035,4.93961};
  
  int n_PbF2=56;
  double en_PbF2[] = {1.55 ,1.569,1.59 ,1.61 ,1.631,1.653,1.675,1.698,1.722,1.746,1.771,1.797,1.823,1.851,1.879,1.907,1.937,1.968,2    ,2.033,2.066,2.101,2.138,2.175,2.214,2.254,2.296,2.339,2.384,2.431,2.48 ,2.53 ,2.583,2.638,2.695,2.755,2.818,2.883,2.952,3.024,3.1  ,3.179,3.263,3.351,3.444,3.542,3.647,3.757,3.875,3.999,4.133,4.275,4.428,4.592,4.769,4.959};
  double ab_PbF2[]= {407,403.3,379.1,406.3,409.7,408.9,406.7,404.7,391.7,397.7,409.6,403.7,403.8,409.7,404.9,404.2,407.1,411.1,403.1,406.1,415.4,399.1,405.8,408.2,385.7,405.6,405.2,401.6,402.6,407.1,417.7,401.1,389.9,411.9,400.9,398.3,402.1,408.7,384.8,415.8,413.1,385.7,353.7,319.1,293.6,261.9,233.6,204.4,178.3,147.6,118.2,78.7 ,51.6 ,41.5 ,24.3 ,8.8};
  double ref_PbF2[]= {1.749,1.749,1.75 ,1.75 ,1.751,1.752,1.752,1.753,1.754,1.754,1.755,1.756,1.757,1.757,1.758,1.759,1.76 ,1.761,1.762,1.764,1.765,1.766,1.768,1.769,1.771,1.772,1.774,1.776,1.778,1.78 ,1.782,1.785,1.787,1.79 ,1.793,1.796,1.8  ,1.804,1.808,1.813,1.818,1.824,1.83 ,1.837,1.845,1.854,1.865,1.877,1.892,1.91 ,1.937,1.991,1.38 ,1.915,1.971,2.019};

  const int n_Sapphire=57;
  double en_Sapphire[] = {1.02212,1.05518,1.09045,1.12610,1.16307,1.20023,1.23984,1.28043,1.32221,1.36561,1.41019,1.45641,1.50393,1.55310,1.60393,1.65643,1.71059,1.76665,1.82437,1.88397,1.94576,2.00946,2.07504,2.14283,2.21321,2.28542,2.36025,2.43727,2.51693,2.59924,2.68480,2.77245,2.86337,2.95693,3.05379,3.15320,3.25674,3.36273,3.47294,3.58646,3.70433,3.82549,3.94979,4.07976,4.21285,4.35032,4.49380,4.64012,4.79258,4.94946,5.11064,5.27816,5.44985,5.62797,5.81266,6.00407,6.19920};
  double ref_Sapphire[] = {1.75188,1.75253,1.75319,1.75382,1.75444,1.75505,1.75567,1.75629,1.75691,1.75754,1.75818,1.75883,1.75949,1.76017,1.76088,1.76160,1.76235,1.76314,1.76395,1.76480,1.76570,1.76664,1.76763,1.76867,1.76978,1.77095,1.77219,1.77350,1.77490,1.77639,1.77799,1.77968,1.78150,1.78343,1.78551,1.78772,1.79011,1.79265,1.79540,1.79835,1.80154,1.80497,1.80864,1.81266,1.81696,1.82163,1.82674,1.83223,1.83825,1.84480,1.85191,1.85975,1.86829,1.87774,1.88822,1.89988,1.91270};

  double ab_Sapphire[n_Sapphire];
  //crystan standard grade 2mm
  double transmitance_Sapphire[]= {0.845,0.844,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.843,0.842,0.842,0.842,0.842,0.838,0.838,0.838,0.838,0.838,0.838,0.836,0.836,0.836,0.836,0.834,0.834,0.833,0.832,0.832,0.831,0.831,0.83,0.828,0.828,0.828,0.828,0.828,0.828,0.828,0.828,0.828,0.825,0.80,0.67,0.47,0.23,0.22};
  for(int i=0;i<n_Sapphire;i++)  ab_Sapphire[i] = (-1)/log(transmitance_Sapphire[i]+2*0.077007)*2*mm;

  
  /*************************** ABSORPTION COEFFICIENTS *****************************/

  // absorption of KamLandOil per 50 cm - from jjv
  double KamLandOilAbsorption[num]=
    {0.97469022,0.976603956,0.978511548,0.980400538,0.982258449,0.984072792,
     0.985831062,0.987520743,0.989129303,0.990644203,0.992052894,
     0.993342822,0.994501428,0.995516151,0.996374433,0.997063719,
     0.997571464,0.997885132,0.997992205,0.997880183,0.997536591,
     0.99,0.98,0.97,0.96,0.94,0.93,0.924507,0.89982,0.883299,
     0.85657,0.842637,0.77020213,0.65727,0.324022,0.019192};

  // absorption of quartz per 1 m - from jjv
  double QuartzAbsorption[num] = 
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541,0.998561611,0.998435332,0.998294892,0.998138345,
     0.997963425,0.997767484,0.997547418,
     0.99729958,0.99701966,0.99670255,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,0.990610945};
  
  // absorption of epotek per one layer - thicknes 0.001'' - from jjv
  double EpotekAbsorption[num] = 
    {0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.9999,0.9998,0.9995,0.999,0.998,0.997,0.996,0.9955,0.993,
     0.9871,0.9745};

  //N-Lak 33a
  double Nlak33aAbsorption[76]={371813,352095,331021,310814,291458,272937,255238,238342,222234,206897,192313,178463,165331,152896,141140,130043,119585,109747,100507,91846.3,83743.1,76176.7,69126.1,62570.2,56488,50858.3,45660.1,40872.4,36474.6,32445.8,28765.9,25414.6,22372.2,19619.3,17136.9,14906.5,12910.2,11130.3,9550.13,8153.3,6924.25,5848.04,4910.46,4098.04,3398.06,2798.54,2288.32,1856.99,1494.92,1193.28,943.973,739.657,573.715,440.228,333.94,250.229,185.064,134.967,96.9664,68.5529,47.6343,32.4882,21.7174,14.2056,9.07612,5.65267,3.4241,2.01226,1.14403,0.62722,0.330414,0.166558,0.0799649,0.0363677,0.0155708,0.00623089};

  double EpotekThickness = 0.001*2.54*cm;
  for(int i=0;i<num;i++){
    WaveLength[i]= (300 +i*10)*nanometer;
    //Absorption[i]= 100*m; // not true, just due to definiton -> not absorb any
    AirAbsorption[i] = 4.*cm; // if photon in the air -> kill it immediately
    AirRefractiveIndex[i] = 1.; 
    PhotonEnergy[num-(i+1)]= LambdaE/WaveLength[i];

    /* as the absorption is given per length and G4 needs 
       mean free path length, calculate it here
       mean free path length - taken as probability equal 1/e
       that the photon will be absorbed */
      
    EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*EpotekThickness;
    QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
    KamLandOilAbsorption[i] = (-1)/log(KamLandOilAbsorption[i])*50*cm;
  }
 

  
  /**************************** REFRACTIVE INDEXES ****************************/
  
  // only phase refractive indexes are necessary -> g4 calculates group itself !!
  
  double QuartzRefractiveIndex[num]={
    1.456535,1.456812,1.4571,1.457399,1.457712,1.458038,1.458378,
    1.458735,1.459108,1.4595,1.459911,1.460344,1.460799,1.46128,
    1.461789,1.462326,1.462897,1.463502,1.464146,1.464833,
    1.465566,1.46635,1.46719,1.468094,1.469066,1.470116,1.471252,1.472485,
    1.473826,1.475289,1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};

  double EpotekRefractiveIndex[num]={
    1.554034,1.555575,1.55698,1.558266,1.559454,1.56056,1.561604,
    1.562604,1.563579,1.564547,1.565526,1.566536,1.567595,
    1.568721,1.569933,1.57125,1.57269,1.574271,1.576012,
    1.577932,1.580049,1.582381,1.584948,1.587768,1.590859,
    1.59424,1.597929,1.601946,1.606307,1.611033,1.616141,1.621651,1.62758,
    1.633947,1.640771,1.64807};

  double KamLandOilRefractiveIndex[num]={
    1.433055,1.433369,1.433698,1.434045,1.434409,1.434793,1.435198,
    1.435626,1.436077,1.436555,1.4371,1.4376,1.4382,1.4388,1.4395,
    1.4402,1.4409,1.4415,1.4425,1.4434,1.4444,1.4455,1.4464,1.4479,1.4501,
    1.450428,1.451976,1.453666,1.455513,1.45754,1.45977,1.462231,1.464958,
    1.467991,1.471377,1.475174};

  double Nlak33aRefractiveIndex[76]={1.73816,1.73836,1.73858,1.73881,1.73904,1.73928,1.73952,1.73976,1.74001,1.74026,1.74052,1.74078,1.74105,1.74132,1.7416,1.74189,1.74218,1.74249,1.74279,1.74311,1.74344,1.74378,1.74412,1.74448,1.74485,1.74522,1.74562,1.74602,1.74644,1.74687,1.74732,1.74779,1.74827,1.74878,1.7493,1.74985,1.75042,1.75101,1.75163,1.75228,1.75296,1.75368,1.75443,1.75521,1.75604,1.75692,1.75784,1.75882,1.75985,1.76095,1.76211,1.76335,1.76467,1.76608,1.76758,1.7692,1.77093,1.77279,1.7748,1.77698,1.77934,1.7819,1.7847,1.78775,1.79111,1.79481,1.79889,1.80343,1.8085,1.81419,1.82061,1.8279,1.83625,1.84589,1.85713,1.87039};

  /* ASSIGNING REFRACTIVE AND ABSORPTION PROPERTIES TO THE GIVEN MATERIALS */

  // Quartz material => Si02
  G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX",       PhotonEnergy, QuartzRefractiveIndex,num);
  QuartzMPT->AddProperty("ABSLENGTH",    PhotonEnergy, QuartzAbsorption,           num);

  // assign this parameter table to BAR material
  BarMaterial->SetMaterialPropertiesTable(QuartzMPT);

  // Air
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX",    PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption,      num);
  //  assign this parameter table to the air 
  defaultMaterial->SetMaterialPropertiesTable(AirMPT);


  // KamLandOil                                                
  G4MaterialPropertiesTable* KamLandOilMPT = new G4MaterialPropertiesTable();
  KamLandOilMPT->AddProperty("RINDEX", PhotonEnergy, KamLandOilRefractiveIndex, num);
  KamLandOilMPT->AddProperty("ABSLENGTH", PhotonEnergy, KamLandOilAbsorption, num);
  // assing this parameter table  to the KamLandOil
  OilMaterial->SetMaterialPropertiesTable(KamLandOilMPT);  

  // N-Lak 33a                                                
  G4MaterialPropertiesTable* Nlak33aMPT = new G4MaterialPropertiesTable();
  Nlak33aMPT->AddProperty("RINDEX", PhotonEnergyNlak33a, Nlak33aRefractiveIndex, 76);
  Nlak33aMPT->AddProperty("ABSLENGTH",PhotonEnergyNlak33a, Nlak33aAbsorption, 76);
  Nlak33aMaterial->SetMaterialPropertiesTable(Nlak33aMPT);

  // PbF2
  G4MaterialPropertiesTable* PbF2MPT = new G4MaterialPropertiesTable();
  PbF2MPT->AddProperty("RINDEX", en_PbF2, ref_PbF2, n_PbF2);
  PbF2MPT->AddProperty("ABSLENGTH",en_PbF2, ab_PbF2, n_PbF2);
  PbF2Material->SetMaterialPropertiesTable(PbF2MPT);

  // Sapphire
  G4MaterialPropertiesTable* SapphireMPT = new G4MaterialPropertiesTable();
  SapphireMPT->AddProperty("RINDEX", en_Sapphire, ref_Sapphire, n_Sapphire);
  SapphireMPT->AddProperty("ABSLENGTH",en_Sapphire, ab_Sapphire, n_Sapphire);
  SapphireMaterial->SetMaterialPropertiesTable(SapphireMPT);

  
  // Epotek Glue                                        
  G4MaterialPropertiesTable* EpotekMPT = new G4MaterialPropertiesTable();
  EpotekMPT->AddProperty("RINDEX", PhotonEnergy, EpotekRefractiveIndex, num);
  EpotekMPT->AddProperty("ABSLENGTH", PhotonEnergy, EpotekAbsorption, num);
  // assign this parameter table to the epotek
  epotekMaterial->SetMaterialPropertiesTable(EpotekMPT);

}


void G4EicDircDetector::SetVisualization(){

  /*G4VisAttributes *waModuleEnvelope = new G4VisAttributes();
  waModuleEnvelope->SetVisibility(false);
  waModuleEnvelope->SetColour(G4Colour::Red());
  waModuleEnvelope->SetForceWireframe(true);
  log_module_envelope->SetVisAttributes(waModuleEnvelope);
  log_module_envelope_inner->SetVisAttributes(waModuleEnvelope);

  G4VisAttributes *waFullEnvelope = new G4VisAttributes();
  waFullEnvelope->SetVisibility(false);
  waFullEnvelope->SetColour(G4Colour::Red());
  waFullEnvelope->SetForceWireframe(true);
  DetectorLog_Det->SetVisAttributes(waFullEnvelope);

  G4VisAttributes *waSupport = new G4VisAttributes();
  waSupport->SetColour(G4Colour::Gray());
  waSupport->SetVisibility(true);
  waSupport->SetForceSolid(true);
  Log_End_Support->SetVisAttributes(waSupport);
  Log_Longitudinal_Support->SetVisAttributes(waSupport);
  Log_End_Support_inner->SetVisAttributes(waSupport);
  */
  G4Colour DircColour = G4Colour(1.,1.0,0.);

  G4VisAttributes *waDirc = new G4VisAttributes(DircColour);
  waDirc->SetForceWireframe(true);
  waDirc->SetVisibility(false);
  lDirc->SetVisAttributes(waDirc);
  
  G4VisAttributes *waFd = new G4VisAttributes(DircColour);
  waFd->SetForceWireframe(true);
  lFd->SetVisAttributes(waFd);

  G4VisAttributes *waBar = new G4VisAttributes(G4Colour(0.,1.,0.9,0.05)); //0.05
  waBar->SetVisibility(true);
  lBarL->SetVisAttributes(waBar);
  lBarS->SetVisAttributes(waBar);
  
  G4VisAttributes *waGlue = new G4VisAttributes(G4Colour(0.,0.4,0.9,0.1));
  waGlue->SetVisibility(true);
  lGlue->SetVisAttributes(waGlue);
  
  G4VisAttributes *waMirror = new G4VisAttributes(G4Colour(1.,1.,0.9,0.2));
  waMirror->SetVisibility(true);
  lMirror->SetVisAttributes(waMirror);

  double transp = 0.4;
  G4VisAttributes * vaLens = new G4VisAttributes(G4Colour(0.,1.,1.,transp));
  vaLens->SetForceWireframe(true);
  // vaLens->SetForceAuxEdgeVisible(true);
  // vaLens->SetForceLineSegmentsPerCircle(10);
  // vaLens->SetLineWidth(4);
   
  if(fLensId==100 )lLens3->SetVisAttributes(vaLens);    
  
  if(fLensId==2 || fLensId==3 || fLensId==6 ){
    lLens1->SetVisAttributes(vaLens);
    G4VisAttributes * vaLens2 = new G4VisAttributes(G4Colour(0.,1.,1.,transp));
    vaLens2->SetColour(G4Colour(0.,0.5,1.,transp));
    vaLens2->SetForceWireframe(true);
    lLens2->SetVisAttributes(vaLens2);
    if(fLensId==3 || fLensId==6) lLens3->SetVisAttributes(vaLens);
  }

  G4VisAttributes *waPrizm = new G4VisAttributes(G4Colour(0.,0.9,0.9,0.4)); //0.4
  // waPrizm->SetForceAuxEdgeVisible(true);
  // waPrizm->SetForceSolid(true);
  lPrizm->SetVisAttributes(waPrizm);

  // G4VisAttributes *waMcp = new G4VisAttributes(G4Colour(0.1,0.1,0.9,0.3));
  G4VisAttributes *waMcp = new G4VisAttributes(G4Colour(1.0,0.,0.1,0.4));
  // waMcp->SetForceWireframe(true);
  waMcp->SetForceSolid(true);
  lMcp->SetVisAttributes(waMcp);

  G4VisAttributes *waPixel = new G4VisAttributes(G4Colour(0.7,0.0,0.1,0.5));
  waPixel->SetForceWireframe(true);
  if(fMcpLayout==3) waPixel->SetVisibility(false);
  else waPixel->SetVisibility(true);
  lPixel->SetVisAttributes(waPixel);

}


void G4EicDircDetector::SetQuantumEfficiency(int id){
  const int num = 36;
  //ideal pmt quantum efficiency
  double QuantumEfficiencyIdial[num]=
    {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0};

  // Burle PMT's 
  double QuantumEfficiencyB[num] =
    {0.,0.001,0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06,
     0.07,0.09,0.1,0.13,0.15,0.17,0.2,0.24,0.26,0.28,0.282,0.284,0.286,
     0.288,0.29,0.28,0.26,0.24,0.22,0.20,0.18,0.15,0.13,0.12,0.10};
  
  //hamamatsu pmt quantum efficiency
  double QuantumEfficiencyPMT[num]=
    {0.001,0.002,0.004,0.007,0.011,0.015,0.020,0.026,0.033,0.040,0.045,
     0.056,0.067,0.085,0.109,0.129,0.138,0.147,0.158,0.170,
     0.181,0.188,0.196,0.203,0.206,0.212,0.218,0.219,0.225,0.230,
     0.228,0.222,0.217,0.210,0.199,0.177};
  
  if(id == 0 ) fQuantumEfficiency = QuantumEfficiencyIdial;
  if(id == 1 ) fQuantumEfficiency = QuantumEfficiencyPMT;
  if(id == 2 ) fQuantumEfficiency = QuantumEfficiencyB;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

}

void G4EicDircDetector::Print(const std::string &what) const
{
  std::cout << "EIC Dirc Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}

