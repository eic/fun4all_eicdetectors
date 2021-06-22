AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4eiccalos.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4eiccalos_la_LIBADD = \
  -lfun4all \
  -lphg4hit \
  -lg4detectors \
  -ltrackbase_historic_io

pkginclude_HEADERS = \
  PHG4BackwardHcalSubsystem.h \
  PHG4CrystalCalorimeterSubsystem.h \
  PHG4ForwardCalCellReco.h \
  PHG4ForwardEcalSubsystem.h \
  PHG4ForwardHcalSubsystem.h \
  PHG4LFHcalSubsystem.h

libg4eiccalos_la_SOURCES = \
  PHG4BackwardHcalDetector.cc \
  PHG4BackwardHcalDisplayAction.cc \
  PHG4BackwardHcalSteppingAction.cc \
  PHG4BackwardHcalSubsystem.cc \
  PHG4CrystalCalorimeterDetector.cc \
  PHG4CrystalCalorimeterDisplayAction.cc \
  PHG4CrystalCalorimeterSteppingAction.cc \
  PHG4CrystalCalorimeterSubsystem.cc \
  PHG4FCalDetector.cc \
  PHG4FCalSteppingAction.cc \
  PHG4FCalSubsystem.cc \
  PHG4ForwardEcalDetector.cc \
  PHG4ForwardEcalDisplayAction.cc \
  PHG4EICForwardEcalDetector.cc \
  PHG4ForwardEcalSteppingAction.cc \
  PHG4ForwardEcalSubsystem.cc \
  PHG4ForwardHcalDetector.cc \
  PHG4ForwardHcalDisplayAction.cc \
  PHG4ForwardHcalSteppingAction.cc \
  PHG4ForwardHcalSubsystem.cc \
  PHG4LFHcalDetector.cc \
  PHG4LFHcalDisplayAction.cc \
  PHG4LFHcalSteppingAction.cc \
  PHG4LFHcalSubsystem.cc \
  PHG4ForwardCalCellReco.cc \
  PHG4FPbScDetector.cc \
  PHG4FPbScSteppingAction.cc \
  PHG4FPbScRegionSteppingAction.cc \
  PHG4FPbScSubsystem.cc \
  PHG4ProjCrystalCalorimeterDetector.cc

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
  libg4eiccalos.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
