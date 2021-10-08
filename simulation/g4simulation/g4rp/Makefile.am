AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libEICG4RP.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -isystem ${G4_MAIN}/include

pkginclude_HEADERS = \
  EICG4RPHitTree.h \
  EICG4RPSubsystem.h

libEICG4RP_la_SOURCES = \
  $(ROOT5_DICTS) \
  EICG4RPSubsystem.cc \
  EICG4RPDetector.cc \
  EICG4RPHitTree.cc \
  EICG4RPSteppingAction.cc

libEICG4RP_la_LDFLAGS  = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lcalo_io \
  -lfun4all \
  -lg4detectors_io \
  -lphg4hit\
  -lgsl \
  -lgslcblas \
  -lg4eval_io \
  -lg4eval \
  -leiceval

libEICG4RP_la_LIBADD = \
  -lphool \
  -lSubsysReco\
  -lg4detectors\
  -lg4testbench 
#   -L/cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/gsl-2.6/lib \ 
#   \ -lm

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD   = libEICG4RP.la


testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
