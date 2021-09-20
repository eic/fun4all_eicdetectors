#ifndef CREATECZHITCONTAINER_H
#define CREATECZHITCONTAINER_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class PHG4TruthInfoContainer;
class PHG4Hit;

class CreateCZHitContainer : public SubsysReco
{
 public:
  //! constructor
  CreateCZHitContainer(const std::string &name = "BMT");

  //! destructor
  virtual ~CreateCZHitContainer();

  //! full initialization
  int InitRun(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  //int End(PHCompositeNode *);
  
  PHG4Hit* merge_hits(PHG4Hit*, PHG4Hit*);

 private:
  std::string _node_postfix;
  PHG4TruthInfoContainer* _truth_container;
  PHG4Hit* _hit_C;
  PHG4Hit* _hit_Z;
  PHG4Hit* _hit_CZ;
};

#endif
