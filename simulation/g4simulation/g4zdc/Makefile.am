AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include\
  -I$(G4_MAIN)/include 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  EICG4ZDCSubsystem.h\
  EICG4ZDCNtuple.h\
  EICG4ZDCHitTree.h\
  EICG4ZDCRawTowerBuilder.h


lib_LTLIBRARIES = \
  libEICG4ZDC.la

libEICG4ZDC_la_SOURCES = \
  EICG4ZDCSubsystem.cc\
  EICG4ZDCDetector.cc\
  EICG4ZDCSteppingAction.cc\
  EICG4ZDCStructure.cc\
  EICG4ZDCNtuple.cc\
  EICG4ZDCHitTree.cc\
  EICG4ZDCRawTowerBuilder.cc

libEICG4ZDC_la_LIBADD = \
  -lphool \
  -lSubsysReco\
  -lg4detectors\
  -lg4testbench\
  -leiczdcbase

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libEICG4ZDC.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
