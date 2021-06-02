#ifndef EICG4ZDCNTUPLE_H
#define EICG4ZDCNTUPLE_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

// Forward declerations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;

class EICG4ZDCNtuple : public SubsysReco
{
 public:
  //! constructor
  EICG4ZDCNtuple(const std::string &name = "EICG4ZDCNtuple", const std::string &filename = "EICG4ZDCNtuple.root");

  //! destructor
  ~EICG4ZDCNtuple() override;

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! end of run method
  int End(PHCompositeNode *) override;

  void AddNode(const std::string &name, const int detid = 0);

 protected:
  int nblocks;
  Fun4AllHistoManager *hm;
  std::string _filename;
  std::set<std::string> _node_postfix;
  std::map<std::string, int> _detid;
  TNtuple *ntup;
  TFile *outfile;
};

#endif
