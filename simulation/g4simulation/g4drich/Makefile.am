AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include \
  -I$(G4_MAIN)/include/Geant4

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  EICG4dRICHSubsystem.h\
  EICG4dRICHTree.h

lib_LTLIBRARIES = \
  libEICG4dRICH.la

libEICG4dRICH_la_SOURCES = \
  EICG4dRICHSubsystem.cc\
  EICG4dRICHDetector.cc\
  EICG4dRICHHit.cc\
  EICG4dRICHSteppingAction.cc\
  EICG4dRICHTree.cc

libEICG4dRICH_la_LIBADD = \
  -lphool \
  -lSubsysReco\
  -lg4detectors\
  -lg4testbench 

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libEICG4dRICH.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
