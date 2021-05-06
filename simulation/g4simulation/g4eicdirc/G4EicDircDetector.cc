#include "G4EicDircDetector.h"

#include "G4EicDircDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
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

#include <cmath>
#include <iostream>  // for operator<<, endl, bas...

class G4VSolid;
class PHCompositeNode;

using namespace std;

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
//_______________________________________________________________
int G4EicDircDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter =
      m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end())
  {
    return 1;
  }

  return 0;
}

void G4EicDircDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  DefineMaterials();
  double fBar[3];
  double  fPrizm[6];
  double fFd[3];
  int fNBar = 1;
  double fMirror[3];
  int  fLensId = 2;
  double fLens[3];
  double fPrizmT[6];
  double fdTilt = 80*deg;
  double fMcpActive[3];
  fPrizmT[0] = 390;
  fPrizmT[1] = 400 - 290*cos(fdTilt); //
  fPrizmT[2] = 290*sin(fdTilt); // hight
  fPrizmT[3] = 50; // face
  fPrizmT[4] = 290;
  fPrizmT[5] = 290*cos(fdTilt);

  fPrizm[0] = 360; 
fPrizm[1] = 300; 
fPrizm[3] = 50; 
fPrizm[2]= fPrizm[3]+300*tan(32*deg);

  fFd[0] = fPrizm[0]; 
fFd[1]=fPrizm[2]; 
fFd[2] =1;

  fBar[0] = 17; 
  fBar[1] = fPrizm[0]/fNBar; 
  fBar[2] = 1050; //4200; //4200

  fMirror[0] = 20; 
fMirror[1] = fPrizm[0]; 
fMirror[2] =1;
fLens[0] = 40.;
fLens[1] = 40; 
fLens[2]=10;
fMcpActive[0] = 53.; 
fMcpActive[1] = 53; 
fMcpActive[2]=1;

  double gluethickness=0.05;
  double dirclength=fBar[2]*4+gluethickness*4;
  G4VSolid* gDirc = new G4Box("gDirc",400.,300.,0.5*dirclength+550);
  G4LogicalVolume *lDirc = new G4LogicalVolume(gDirc,G4Material::GetMaterial("G4_AIR"),"lDirc",0,0,0);
  G4VSolid* gFd = new G4Box("gFd",0.5*fFd[1],0.5*fFd[0],0.5*fFd[2]);
   G4LogicalVolume *lFd = new G4LogicalVolume(gFd,G4Material::GetMaterial("G4_AIR"),"lFd",0,0,0);
   int fNBoxes = 1;
  double fRadius = 970;
  double tphi, dphi = 22.5*deg; //22.5*deg;
  int BoxId = 0;
 G4PVPlacement *phy = new G4PVPlacement(0,G4ThreeVector(0,0,0),lDirc,"wDirc",logicWorld,false,0);
  m_PhysicalVolumesSet.insert(phy);
  // The Bar
  G4VSolid *gBar = new G4Box("gBar",fBar[0]/2.,fBar[1]/2.,fBar[2]/2.);
  G4LogicalVolume *lBar = new G4LogicalVolume(gBar,G4Material::GetMaterial("quartz"),"lBar",0,0,0);
  G4VSolid *gExpVol = new G4Box("gExpVol",fBar[0]/2.,fPrizm[0]/2.,fBar[2]/2.);
   G4LogicalVolume *lExpVol = new G4LogicalVolume(gExpVol,G4Material::GetMaterial("quartz"),"lExpVol",0,0,0);
  G4VSolid *gGlue = new G4Box("gGlue",fBar[0]/2.,fBar[1]/2.,0.5*gluethickness);
  G4LogicalVolume *lGlue = new G4LogicalVolume(gGlue,G4Material::GetMaterial("Epotek"),"lGlue",0,0,0);
  G4VSolid *gGlueE = new G4Box("gGlueE",fBar[0]/2.,fPrizm[0]/2.,0.5*gluethickness);
  G4LogicalVolume *lGlueE = new G4LogicalVolume(gGlueE,G4Material::GetMaterial("Epotek"),"lGlueE",0,0,0);
  if(fNBar==1){
    for(int j=0;j<4; j++){
      double z = -0.5*dirclength+0.5*1050+(1050+gluethickness)*j;
      phy =  new G4PVPlacement(0,G4ThreeVector(0,0,z),lBar,"wBar", lDirc,false,0);
  m_PhysicalVolumesSet.insert(phy);

      phy =  new G4PVPlacement(0,G4ThreeVector(0,0,z+0.5*(1050+gluethickness)),lGlue,"wGlue", lDirc,false,0);
  m_PhysicalVolumesSet.insert(phy);
    }
  }
  // The Mirror
  G4Box* gMirror = new G4Box("gMirror",fMirror[0]/2.,fMirror[1]/2.,fMirror[2]/2.);
  G4LogicalVolume * lMirror = new G4LogicalVolume(gMirror,G4Material::GetMaterial("G4_Al"),"lMirror",0,0,0);
  phy =new G4PVPlacement(0,G4ThreeVector(0,0,-0.5*dirclength-fMirror[2]/2.),lMirror,"wMirror", lDirc,false,0);
  m_PhysicalVolumesSet.insert(phy);

  // The Lenses
  if(fLensId == 2){ // 2-layer lens
    G4double lensrad = 70;
    G4double lensMinThikness = 2;
    G4VSolid *gfbox = new G4Box("Fbox",fLens[0]/2.,fLens[1]/2.,fLens[2]/2.);
    G4ThreeVector zTrans(0, 0, -lensrad+fLens[2]/2.-lensMinThikness);

    G4Sphere* gsphere = new G4Sphere("Sphere",0,70,0.*deg,360.*deg,0.*deg,380.*deg);
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Fbox*Sphere", gfbox, gsphere,new G4RotationMatrix(),zTrans); 
    G4SubtractionSolid* gLens2 = new G4SubtractionSolid("Fbox-Sphere", gfbox, gsphere,new G4RotationMatrix(),zTrans);

    G4LogicalVolume *lLens1 = new G4LogicalVolume(gLens1,G4Material::GetMaterial("Nlak33a"),"lLens1",0,0,0); //Nlak33aMaterial  
     G4LogicalVolume *lLens2 = new G4LogicalVolume(gLens2,G4Material::GetMaterial("quartz"),"lLens2",0,0,0);
  }
  // The Prizm
  G4VSolid* gPrizm = new G4Trap("gPrizm",fPrizm[0],fPrizm[1],fPrizm[2],fPrizm[3]);
  G4LogicalVolume *lPrizm = new G4LogicalVolume(gPrizm, G4Material::GetMaterial("quartz"),"lPrizm",0,0,0);
  G4VSolid* gPrizmT1 = new G4Trap("gPrizmT1",fPrizmT[0],fPrizmT[1],fPrizmT[2],fPrizmT[3]);

  G4LogicalVolume *lPrizmT1 = new G4LogicalVolume(gPrizmT1, G4Material::GetMaterial("quartz"),"lPrizmT1",0,0,0);
  
  G4VSolid* gPrizmT2 = new G4Trap("gPrizmT2",fPrizmT[0],fPrizmT[5],fPrizmT[2],0.0001);
  G4LogicalVolume *lPrizmT2 = new G4LogicalVolume(gPrizmT2, G4Material::GetMaterial("quartz"),"lPrizmT2",0,0,0);

  G4RotationMatrix* xRot = new G4RotationMatrix();
  xRot->rotateX(M_PI/2.*rad);

  G4RotationMatrix* fdRot = new G4RotationMatrix();
  G4RotationMatrix *fdrot = new G4RotationMatrix();
  G4double evshiftz = 0.5*dirclength+fPrizm[1]+fMcpActive[2]/2.+fLens[2];
  G4double evshiftx = 0;
  return;
}

void G4EicDircDetector::DefineMaterials()
{
  double density;
  double a;
  double z;
  int ncomponents;
  int natoms;
 G4String symbol;
  G4Element* Si = new G4Element("Silicon" ,symbol="Si", z= 14., a= 28.09*g/mole);
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  // quartz material = SiO2
  G4Material* SiO2 = new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);
  G4Material* Epotek = new G4Material("Epotek",density=1.2*g/cm3,ncomponents=3);

  Epotek->AddElement(C,natoms=3);
  Epotek->AddElement(H,natoms=5);
  Epotek->AddElement(O,natoms=2);

  G4Material* Nlak33aMaterial  = new G4Material("Nlak33a",density= 4.220*g/cm3, ncomponents=2);
  Nlak33aMaterial->AddElement(Si, natoms=1);
  Nlak33aMaterial->AddElement(O , natoms=2);


}
void G4EicDircDetector::Print(const std::string &what) const
{
  cout << "EIC Dirc Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Version 0.1" << endl;
    cout << "Parameters:" << endl;
    m_Params->Print();
  }
  return;
}
