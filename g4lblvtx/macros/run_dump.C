#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <nodedump/Dumper.h>
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libphnodedump.so)
#endif
void run_dump(const char *infile, const int evts=100)
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("libphnodedump.so");
  gSystem->Load("libg4dst.so");

  Fun4AllServer* se = Fun4AllServer::instance();

  Dumper *dmp = new Dumper();
  gSystem->Exec("mkdir ./asciidump");
  dmp->SetOutDir("./asciidump");

  se->registerSubsystem(dmp);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  se->registerInputManager(in);
  se->fileopen("DSTin",infile);
  se->run(evts);
  se->End();
  delete se;
}
