#ifndef G4DIRCTREE_H
#define G4DIRCTREE_H

#include "G4EventTree.h"

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <TVector3.h>

// Forward declarations
class PHCompositeNode;
class PHG4HitContainer;
class TFile;
class TTree;

class G4DIRCTree : public SubsysReco
{
 public:
  //! constructor
  G4DIRCTree(const std::string &name = "G4DIRCTree", const std::string &filename = "G4DIRCTree.root");

  //! destructor
  ~G4DIRCTree() override {}

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! hit processing method
  int process_hit(PHG4HitContainer *hits, const std::string &dName, int detid, int &nhits);

  //! end of run method
  int End(PHCompositeNode *) override;

  void AddNode(const std::string &name, const int detid = 0);

 protected:
  int nblocks;
  //  std::vector<TH2 *> nhit_edep;
  std::string _filename;
  std::set<std::string> _node_postfix;
  std::map<std::string, int> _detid;

  TTree *g4tree;
  G4EventTree mG4EvtTree;
  TFile *outfile;
};

#endif
