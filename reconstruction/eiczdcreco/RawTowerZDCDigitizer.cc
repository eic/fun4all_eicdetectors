#include "RawTowerZDCDigitizer.h"

#include <eiczdcbase/RawTowerZDC.h>  // for RawTowerZDC
#include <eiczdcbase/RawTowerZDCContainer.h>
#include <eiczdcbase/RawTowerZDCDeadMap.h>
#include <eiczdcbase/RawTowerZDCDefs.h>  // for keytype
#include <eiczdcbase/RawTowerZDCGeom.h>
#include <eiczdcbase/RawTowerZDCGeomContainer.h>
#include <eiczdcbase/RawTowerZDCv1.h>

#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBOSITY_MORE
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phparameter/PHParameters.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cmath>
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>  // for pair

using namespace std;

RawTowerZDCDigitizer::RawTowerZDCDigitizer(const std::string &name)
  : SubsysReco(name)
  , m_DigiAlgorithm(kNo_digitization)
  , m_SimTowers(nullptr)
  , m_RawTowers(nullptr)
  , m_RawTowerGeom(nullptr)
  , m_DeadMap(nullptr)
  , m_Detector("NONE")
  , m_SimTowerNodePrefix("SIM")
  , m_RawTowerNodePrefix("RAW")
  , m_PhotonElecYieldVisibleGeV(NAN)
  , m_PhotonElecADC(NAN)
  , m_PedestalCentralADC(NAN)
  , m_PedestalWidthADC(NAN)
  , m_pedestalFile(false)
  , m_ZeroSuppressionADC(-1000)  //default to apply no zero suppression
  , m_ZeroSuppressionFile(false)
  , m_TowerType(-1)
  , m_SiPMEffectivePixel(40000 * 4)  // sPHENIX EMCal default, 4x Hamamatsu S12572-015P MPPC [sPHENIX TDR]
  , _tower_params(name)
{
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  m_Seed = PHRandomSeed();  // fixed seed handled in PHRandomSeed()
  // cout << Name() << " Random Seed: " << m_Seed << endl;
  gsl_rng_set(m_RandomGenerator, m_Seed);
}

RawTowerZDCDigitizer::~RawTowerZDCDigitizer()
{
  gsl_rng_free(m_RandomGenerator);
}

void RawTowerZDCDigitizer::set_seed(const unsigned int iseed)
{
  m_Seed = iseed;
  gsl_rng_set(m_RandomGenerator, m_Seed);
}

int RawTowerZDCDigitizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "DST Node missing, doing nothing." << endl;
    exit(1);
  }

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    cout << e.what() << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerZDCDigitizer::process_event(PHCompositeNode */*topNode*/)
{
  if (Verbosity())
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "Process event entered. "
         << "Digitalization method: ";

    if (m_DigiAlgorithm == kNo_digitization)
    {
      cout << "directly pass the energy of sim tower to digitalized tower";
    }
    else if (m_DigiAlgorithm == kSimple_photon_digitization)
    {
      cout << "simple digitization with photon statistics, ADC conversion and pedestal";
    }
    cout << endl;
  }
  // loop over all possible towers, even empty ones. The digitization can add towers containing
  // pedestals
  RawTowerZDCGeomContainer::ConstRange all_towers = m_RawTowerGeom->get_tower_geometries();

  double deadChanEnergy = 0;

  for (RawTowerZDCGeomContainer::ConstIterator it = all_towers.first;
       it != all_towers.second; ++it)
  {
    const RawTowerZDCDefs::keytype key = it->second->get_id();
    //    RawTowerZDCDefs::CalorimeterId caloid = RawTowerZDCDefs::decode_caloid(key);
    const int eta = it->second->get_bineta();
    const int phi = it->second->get_binphi();
    const int twr = it->second->get_binl();

    if (m_ZeroSuppressionFile == true)
    {
      const string zsName = "ZS_ADC_eta" + to_string(eta) + "_phi" + to_string(phi) + "_twr" + to_string(twr);
      m_ZeroSuppressionADC = _tower_params.get_double_param(zsName);
    }

    if (m_pedestalFile == true)
    {
      const string pedCentralName = "PedCentral_ADC_eta" + to_string(eta) + "_phi" + to_string(phi) + "_twr" + to_string(twr);
      m_PedestalCentralADC = _tower_params.get_double_param(pedCentralName);
      const string pedWidthName = "PedWidth_ADC_eta" + to_string(eta) + "_phi" + to_string(phi) + "_twr" + to_string(twr);
      m_PedestalWidthADC = _tower_params.get_double_param(pedWidthName);
    }

    if (m_TowerType >= 0)
    {
      // Skip towers that don't match the type we are supposed to digitize
      if (m_TowerType != it->second->get_tower_type())
      {
        continue;
      }
    }

    RawTowerZDC *sim_tower = m_SimTowers->getTower(key);
    if (m_DeadMap)
    {
      if (m_DeadMap->isDeadTower(key))
      {
        if (sim_tower) deadChanEnergy += sim_tower->get_energy();

        sim_tower = nullptr;

        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
               << " apply dead tower " << key << endl;
        }
      }
    }

    RawTowerZDC *digi_tower = nullptr;

    if (m_DigiAlgorithm == kNo_digitization)
    {
      // for no digitization just copy existing towers
      if (sim_tower)
      {
        digi_tower = new RawTowerZDCv1(*sim_tower);
      }
    }
    else if (m_DigiAlgorithm == kSimple_photon_digitization)
    {
      // for photon digitization towers can be created if sim_tower is null pointer
      digi_tower = simple_photon_digitization(sim_tower);
    }
    else if (m_DigiAlgorithm == kSiPM_photon_digitization)
    {
      // for photon digitization towers can be created if sim_tower is null pointer
      digi_tower = sipm_photon_digitization(sim_tower);
    }
    else
    {
      cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
           << " invalid digitization algorithm #" << m_DigiAlgorithm
           << endl;

      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (digi_tower)
    {
      m_RawTowers->AddTower(key, digi_tower);

      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
             << " output tower:"
             << endl;
        digi_tower->identify();
      }
    }
  }

  if (Verbosity())
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "input sum energy = " << m_SimTowers->getTotalEdep() << " GeV"
         << ", dead channel masked energy = " << deadChanEnergy << " GeV"
         << ", output sum digitalized value = " << m_RawTowers->getTotalEdep() << " ADC"
         << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

RawTowerZDC *
RawTowerZDCDigitizer::simple_photon_digitization(RawTowerZDC *sim_tower)
{
  RawTowerZDC *digi_tower = nullptr;

  double energy = 0;
  if (sim_tower)
  {
    energy = sim_tower->get_energy();
  }
  const double photon_count_mean = energy * m_PhotonElecYieldVisibleGeV;
  const int photon_count = gsl_ran_poisson(m_RandomGenerator, photon_count_mean);
  const int signal_ADC = floor(photon_count / m_PhotonElecADC);

  const double pedestal = m_PedestalCentralADC + ((m_PedestalWidthADC > 0) ? gsl_ran_gaussian(m_RandomGenerator, m_PedestalWidthADC) : 0);
  const int sum_ADC = signal_ADC + (int) pedestal;

  if (sum_ADC > m_ZeroSuppressionADC)
  {
    // create new digitalizaed tower
    if (sim_tower)
    {
      digi_tower = new RawTowerZDCv1(*sim_tower);
    }
    else
    {
      digi_tower = new RawTowerZDCv1();
    }
    digi_tower->set_energy((double) sum_ADC);
  }

  if (Verbosity() >= 2)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << endl;

    cout << "input: ";
    if (sim_tower)
    {
      sim_tower->identify();
    }
    else
    {
      cout << "None" << endl;
    }
    cout << "output based on "
         << "sum_ADC = " << sum_ADC << ", zero_sup = "
         << m_ZeroSuppressionADC << " : ";
    if (digi_tower)
    {
      digi_tower->identify();
    }
    else
    {
      cout << "None" << endl;
    }
  }

  return digi_tower;
}

RawTowerZDC *
RawTowerZDCDigitizer::sipm_photon_digitization(RawTowerZDC *sim_tower)
{
  RawTowerZDC *digi_tower = nullptr;

  double energy = 0;
  if (sim_tower)
  {
    energy = sim_tower->get_energy();
  }

  int signal_ADC = 0;

  if (energy > 0)
  {
    const double photon_count_mean = energy * m_PhotonElecYieldVisibleGeV;
    const double poission_param_per_pixel = photon_count_mean / m_SiPMEffectivePixel;
    const double prob_activated_per_pixel = gsl_cdf_poisson_Q(0, poission_param_per_pixel);
    const double active_pixel = gsl_ran_binomial(m_RandomGenerator, prob_activated_per_pixel, m_SiPMEffectivePixel);
    signal_ADC = floor(active_pixel / m_PhotonElecADC);
  }

  const double pedestal = m_PedestalCentralADC + ((m_PedestalWidthADC > 0) ? gsl_ran_gaussian(m_RandomGenerator, m_PedestalWidthADC) : 0);
  const int sum_ADC = signal_ADC + (int) pedestal;

  if (sum_ADC > m_ZeroSuppressionADC)
  {
    // create new digitalizaed tower
    if (sim_tower)
    {
      digi_tower = new RawTowerZDCv1(*sim_tower);
    }
    else
    {
      digi_tower = new RawTowerZDCv1();
    }
    digi_tower->set_energy((double) sum_ADC);
  }

  if (Verbosity() >= 2)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << endl;

    cout << "input: ";
    if (sim_tower)
    {
      sim_tower->identify();
    }
    else
    {
      cout << "None" << endl;
    }
    cout << "output based on "
         << "sum_ADC = " << sum_ADC << ", zero_sup = "
         << m_ZeroSuppressionADC << " : ";
    if (digi_tower)
    {
      digi_tower->identify();
    }
    else
    {
      cout << "None" << endl;
    }
  }

  return digi_tower;
}

void RawTowerZDCDigitizer::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "Run Node missing, doing nothing." << endl;
    throw std::runtime_error(
        "Failed to find Run node in RawTowerZDCDigitizer::CreateNodes");
  }

  string TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  m_RawTowerGeom = findNode::getClass<RawTowerZDCGeomContainer>(topNode, TowerGeomNodeName);
  if (!m_RawTowerGeom)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << " " << TowerGeomNodeName << " Node missing, doing bail out!"
         << endl;
    throw std::runtime_error("Failed to find " + TowerGeomNodeName + " node in RawTowerZDCDigitizer::CreateNodes");
  }

  if (Verbosity() >= 1)
  {
    m_RawTowerGeom->identify();
  }

  const string deadMapName = "DEADMAP_" + m_Detector;
  m_DeadMap = findNode::getClass<RawTowerZDCDeadMap>(topNode, deadMapName);
  if (m_DeadMap)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << " use dead map: ";
    m_DeadMap->identify();
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "DST Node missing, doing nothing." << endl;
    throw std::runtime_error("Failed to find DST node in RawTowerZDCDigitizer::CreateNodes");
  }

  string SimTowerNodeName = "TOWER_" + m_SimTowerNodePrefix + "_" + m_Detector;
  m_SimTowers = findNode::getClass<RawTowerZDCContainer>(dstNode, SimTowerNodeName);
  if (!m_SimTowers)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << " " << SimTowerNodeName << " Node missing, doing bail out!"
         << endl;
    throw std::runtime_error("Failed to find " + SimTowerNodeName + " node in RawTowerZDCDigitizer::CreateNodes");
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }

  // Be careful as a previous digitizer may have been registered for this detector
  string RawTowerNodeName = "TOWER_" + m_RawTowerNodePrefix + "_" + m_Detector;
  m_RawTowers = findNode::getClass<RawTowerZDCContainer>(DetNode, RawTowerNodeName);
  if (!m_RawTowers)
  {
    m_RawTowers = new RawTowerZDCContainer(m_SimTowers->getCalorimeterID());
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_RawTowers,
                                                                   RawTowerNodeName, "PHObject");
    DetNode->addNode(towerNode);
  }

  return;
}
