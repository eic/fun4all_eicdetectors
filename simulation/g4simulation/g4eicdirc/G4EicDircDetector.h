// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4EICDIRCDETECTOR_H
#define G4EICDIRCDETECTOR_H

#include <g4main/PHG4Detector.h>
#include <Geant4/G4Types.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4LogicalVolume.hh>

#include <set>
#include <map>
#include <string>  // for string

class G4EicDircDisplayAction;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class G4EicDircDetector : public PHG4Detector
{
 public:
  //! constructor
  G4EicDircDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~G4EicDircDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  void SetVisualization();
  void SetQuantumEfficiency(G4int id);

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  //@}
 
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  std::string name_base = "test";

  void
    OverlapCheck(bool check = true) override
  {
    overlapcheck_sector = check;
  }

 
 private:
  G4LogicalVolume* lFront;
  G4LogicalVolume* lDirc;
  G4LogicalVolume* lFd;
  G4LogicalVolume* lBar;
  G4LogicalVolume* lGlue;
  G4LogicalVolume* lMirror;
  G4LogicalVolume* lLens1;
  G4LogicalVolume* lLens2;
  G4LogicalVolume* lLens3;
  G4LogicalVolume* lPrizm;
  G4LogicalVolume* lMcp;
  G4LogicalVolume* lPixel;
  G4LogicalVolume* lExpVol;
  G4LogicalVolume* lGlueE;
  G4VPhysicalVolume*   pDirc[100];

  G4VPhysicalVolume* wBar;
  G4VPhysicalVolume* wGlue;
  G4VPhysicalVolume* wMirror;
  G4VPhysicalVolume* wDirc;

  G4Material*        defaultMaterial; // material for bars
  G4Material*        BarMaterial; // material for bars
  G4Material*        OilMaterial;
  G4Material*        MirrorMaterial; // material of mirror
  G4Material*        epotekMaterial;  
  G4Material*        Nlak33aMaterial;
  G4Material*        PbF2Material;
  G4Material*        SapphireMaterial;
  G4Material*        frontMaterial;
  
  G4int fNRow;
  G4int fNCol;
  G4int fNBoxes;
  G4double fRadius;
  G4double fNpix1;
  G4double fNpix2;
  G4double fBoxWidth;
  G4int fGeomType;
  G4int fMcpLayout;
  G4int fLensId;
  G4double fNBar;
  G4double fBar[3];
  G4double fMirror[3];
  G4double fFd[3];
  G4double fPrizm[4];
  G4double fLens[4];
  G4double fMcpTotal[3];
  G4double fMcpActive[3];
  G4ThreeVector fPrismShift;
  G4double fBarsGap;

  G4double fRotAngle;
  G4RotationMatrix *fPrtRot;
  G4double *fQuantumEfficiency;


 protected:
  void DefineMaterials();
  PHParameters *m_Params;

  G4EicDircDisplayAction *m_DisplayAction;
  bool overlapcheck_sector;
  
  // active volumes
  G4LogicalVolume* RegisterLogicalVolume(G4LogicalVolume *);
  G4PVPlacement* RegisterPhysicalVolume(G4PVPlacement *v, const bool active = false);

  std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;
  std::map<G4VPhysicalVolume *, int> m_PhysicalVolumes_active;
  //std::map<G4LogicalVolume*, int> m_logVol_active;
  std::string m_SuperDetector;
};

#endif  // G4EICDIRCDETECTOR_H
