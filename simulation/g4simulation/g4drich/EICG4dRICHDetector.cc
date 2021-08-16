#include "EICG4dRICHDetector.h"
#include "EICG4dRICHConfig.hh"
#include "EICG4dRICHOptics.hh"

#include <fun4all/Fun4AllBase.h>

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>

#include <G4Box.hh>
#include <G4Color.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>
#include <G4tgbVolumeMgr.hh>

#include <cmath>
#include <iostream>

class G4VSolid;
class PHCompositeNode;

// ---------------------------------------------------
EICG4dRICHDetector::EICG4dRICHDetector(PHG4Subsystem *subsys, PHCompositeNode *Node,
                                       PHParameters *parameters, const std::string &dnam)
    : PHG4Detector(subsys, Node, dnam)
    , m_Params(parameters)
    {}

// ---------------------------------------------------
int EICG4dRICHDetector::IsInDetector(G4VPhysicalVolume *volume) const 
{
  std::set<G4VPhysicalVolume *>::const_iterator iter =
      m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end()) 
  {
    return 1;
  }
  return 0;
}

// ---------------------------------------------------
void EICG4dRICHDetector::ConstructMe(G4LogicalVolume *logicWorld) 
{
  EICG4dRICHConfig cfg;
  cfg.model_file = m_Params->get_string_param("mapping_file");
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[+] MODEL TEXT FILE: " << cfg.model_file << std::endl;
  // - check existence
  std::ifstream mf(cfg.model_file.data());
  if (!mf.is_open()) 
  {
    std::cerr << "[+] ERROR in " << __FILE__ << ": cannot find MODEL TEXT FILE" << std::endl;
    return;
  }
  mf.close();

  // graphics
  G4VisAttributes *vis = new G4VisAttributes(G4Color(0., 0., 0.9, 0.5));
  vis->SetForceSolid(true);
  vis->SetVisibility(true);
  vis->SetLineWidth(1);
  vis->SetForceSolid(true);
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  // logicWorld->SetVisAttributes(vis);

  // build detector by text file
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[+] read model text file" << std::endl;
  G4tgbVolumeMgr *volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(cfg.model_file);
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) std::cout << "[+] construct detector from text file" << std::endl;
  G4VPhysicalVolume *vesselPhysVol = volmgr->ReadAndConstructDetector();
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE) 
  {
    std::cout << "[+] detector summary" << std::endl;
    volmgr->DumpSummary();
    std::cout << "[+] detector G4Solid list" << std::endl;
    volmgr->DumpG4SolidList();
  }

  // material optical properties (see shared header EICG4dRICHOptics.hh)
  // - aerogel
  auto aeroPO = new EICG4dRICHAerogel("EICG4dRICHaerogelMat");
  aeroPO->setOpticalParams(cfg.aerOptModel); // mode=3: use experimental data
  // - acrylic filter
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
  {
    std::cout << "[+] Acrylic Wavelength Threshold : " << cfg.filter_thr / (nm) << " nm\n" << std::endl;
  }
  auto acryPO = new EICG4dRICHFilter("EICG4dRICHfilterMat");
  acryPO->setOpticalParams(cfg.filter_thr);
  // - gas radiator
  auto gasPO = new EICG4dRICHGas("EICG4dRICHgasMat");
  gasPO->setOpticalParams();
  // - photo sensors
  auto photoSensor = new EICG4dRICHPhotosensor("EICG4dRICHpsst");
  photoSensor->setOpticalParams("EICG4dRICH");
  // - mirror (simular to photosensor, but different params)
  auto mirror = new EICG4dRICHMirror("EICG4dRICHmirror");
  mirror->setOpticalParams("EICG4dRICH");

  // add to logical world
  logicWorld->AddDaughter(vesselPhysVol);

  // activate volumes, for hit readout
  this->ActivateVolumeTree(vesselPhysVol);

  return;
}

// ---------------------------------------------------
// recursively add detectors to active volume list, descending the tree
// from `volu`
// - use the "activation filter" to decide for which volumes to save hits
// - the petal number is added to `m_PetalMap`, which together
//   with the copy number, provides a unique ID for each photo sensor
void EICG4dRICHDetector::ActivateVolumeTree(G4VPhysicalVolume *volu, G4int petal) 
{

  // get objects
  G4String voluName = volu->GetName();
  G4int voluCopyNo = volu->GetCopyNo();
  G4LogicalVolume *logi = volu->GetLogicalVolume();

  // obtain petal number
  if (voluName.contains("petal"))
  {
    petal = voluCopyNo;
  }

  // activation filter: use this to decide which volumes to save
  // hits for, i.e., which volumes are "active"
  // TODO: need to decide what volume we want to be active
  G4bool activate = true;// voluName.contains("psst") || voluName.contains("vessel");
  if (activate) 
  {
    if (Verbosity() >= Fun4AllBase::VERBOSITY_SOME)  
    {
      std::cout << "[+] activate " << voluName << " petal " << petal << " copy "
           << voluCopyNo << std::endl;
    }
    m_PhysicalVolumesSet.insert(volu);
    m_PetalMap.insert(std::pair<G4VPhysicalVolume *, G4int>(volu, petal));
  }

  // loop over daughters
  G4int nd = logi->GetNoDaughters();
  for (int d = 0; d < nd; d++) 
  {
    this->ActivateVolumeTree(logi->GetDaughter(d), petal);
  }
}

// ---------------------------------------------------
// get petal number
int EICG4dRICHDetector::GetPetal(G4VPhysicalVolume *volu) 
{
  int petalNum;
  try 
  {
    petalNum = m_PetalMap.at(volu);
  } 
  catch (const std::out_of_range &ex) 
  {
    std::cerr << "ERROR in EICG4dRICHDetector: cannot find petal associated with volume"
         << std::endl;
    return -1;
  }
  return petalNum;
}

// ---------------------------------------------------
// get PSST number
int EICG4dRICHDetector::GetPSST(G4VPhysicalVolume *volu) 
{
  return volu->GetName().contains("psst") ? volu->GetCopyNo() : 0;
}

// ---------------------------------------------------
void EICG4dRICHDetector::Print(const std::string &what) const 
{
  std::cout << "EICG4dRICH Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME") 
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}
