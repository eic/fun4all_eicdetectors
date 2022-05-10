AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libbst.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libbst_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lphg4hit \
  -lg4testbench \
  -lg4detectors \
  -lg4detectors_io \
  -ltrackbase_historic_io \
  -lphparameter \
  -lgsl \
  -lgslcblas \
  -lcalo_io

pkginclude_HEADERS = \
  PHG4BSTSubsystem.h

libbst_la_SOURCES = \
  PHG4BSTDetector.cc \
  PHG4BSTDisplayAction.cc \
  PHG4BSTSteppingAction.cc \
  PHG4BSTSubsystem.cc

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
  libbst.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
