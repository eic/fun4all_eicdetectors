AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4eicdirc.la

AM_CPPFLAGS = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I$(G4_MAIN)/include \
    -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4eicdirc_la_LIBADD = \
  -lSubsysReco \
  -lg4detectors \
  -lg4testbench 

pkginclude_HEADERS = \
  G4EicDircSubsystem.h \
  PrtOpBoundaryProcess.h \
  PrtHit.h \
  G4EventTree.h \
  G4DIRCTree.h 

ROOT_IO_DICTS = \
  PrtHit_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  PrtHit_Dict_rdict.pcm 

libg4eicdirc_la_SOURCES = \
  G4EicDircDetector.cc \
  G4EicDircDisplayAction.cc \
  G4EicDircOpBoundaryProcess.cc \
  G4EicDircSubsystem.cc \
  G4EicDircStackingAction.cc \
  G4EicDircSteppingAction.cc \
  PrtHit.cxx \
  PrtHit_Dict.cc \
  PrtOpBoundaryProcess.cxx \
  G4DIRCTree.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4eicdirc


testexternals_g4eicdirc_SOURCES = testexternals.cc
testexternals_g4eicdirc_LDADD = libg4eicdirc.la


testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
