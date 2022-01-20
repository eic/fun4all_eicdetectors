##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libeiczdcreco.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I${G4_MAIN}/include

pkginclude_HEADERS = \
  RawTowerZDCDigitizer.h \
  RawTowerZDCCalibration.h
libeiczdcreco_la_SOURCES = \
  RawTowerZDCDigitizer.cc \
  RawTowerZDCCalibration.cc

libeiczdcreco_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libeiczdcreco_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lg4detectors_io \
  -lg4testbench \
  -lphg4hit \
  -lphparameter \
  -lgsl \
  -lgslcblas \
  -leiczdcbase

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_eiczdcreco

BUILT_SOURCES = testexternals.cc

testexternals_eiczdcreco_SOURCES = testexternals.cc
testexternals_eiczdcreco_LDADD = libeiczdcreco.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)
