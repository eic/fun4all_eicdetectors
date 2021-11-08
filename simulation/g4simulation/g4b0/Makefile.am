AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libEICG4B0.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(ROOTSYS)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/include/eigen3  \
  -I$(ROOTSYS)/include \
  -I${G4_MAIN}/include

pkginclude_HEADERS = \
  EICG4B0Subsystem.h \
  B0TrackFastSimEval.h\
  B0TrackFastSim.h

libEICG4B0_la_SOURCES = \
  $(ROOT5_DICTS) \
  EICG4B0Subsystem.cc \
  EICG4B0Detector.cc \
  EICG4B0SteppingAction.cc \
  B0TrackFastSimEval.cc \ 
  B0TrackFastSim.cc 

libEICG4B0_la_LDFLAGS  = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
libEICG4B0_la_LIBADD = \
  -lcalo_io \
  -lfun4all \
  -lPHGenFit \
  -lg4detectors_io \
  -lphg4hit\
  -lgsl \
  -lgslcblas \
  -lg4eval_io \
  -lg4eval \
  -leiceval\
  -lgslcblas \
  -lphgeom \
  -lphg4hit \
  -ltrackbase_historic_io \
  -lphparameter\
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
testexternals_LDADD   = libEICG4B0.la


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
