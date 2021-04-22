// assuming the file contains our PHGeometry
void TGeoToGdml(const char *infile="genfitGeom_AllSi_v3.root", const char *tgeoobj="FAIRGeom")
{
  TFile *f = TFile::Open(infile);
  TGeoManager *geomanager = ( TGeoManager *) f->Get(tgeoobj);
  char outfilename[200];
  sprintf(outfilename,"%s.gdml",tgeoobj);
  geomanager->Export(outfilename);
  gSystem->Exit(0);
}
