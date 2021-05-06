AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4lblvtx.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4lblvtx_la_LIBADD = \
  -lfun4all \
  -lphg4hit \
  -lg4detectors \
  -ltrackbase_historic_io

pkginclude_HEADERS = \
  AllSiliconTrackerSubsystem.h \
  EicFRichSubsystem.h \
  AllSi_Al_support_Subsystem.h \
  G4LBLVtxSubsystem.h \
  SimpleNtuple.h \
  TrackFastSimEval.h \
  PHG4ParticleGenerator_flat_pT.h

libg4lblvtx_la_SOURCES = \
  AllSiliconTrackerDetector.cc \
  AllSiliconTrackerDisplayAction.cc \
  AllSiliconTrackerSteppingAction.cc \
  AllSiliconTrackerSubsystem.cc \
  EicFRichDetector.cc \
  EicFRichSteppingAction.cc \
  EicFRichSubsystem.cc \
  AllSi_Al_support_Detector.cc \
  AllSi_Al_support_SteppingAction.cc \
  AllSi_Al_support_Subsystem.cc \
  AllSi_Al_support_DisplayAction.cc \
  G4LBLVtxDetector.cc \
  G4LBLVtxDisplayAction.cc \
  G4LBLVtxSteppingAction.cc \
  G4LBLVtxSubsystem.cc \
  SimpleNtuple.cc \
  TrackFastSimEval.cc \
  PHG4ParticleGenerator_flat_pT.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libg4lblvtx.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
