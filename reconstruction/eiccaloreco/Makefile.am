AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  RawClusterBuilderkV3.h

lib_LTLIBRARIES = \
  libEICCaloReco.la

libEICCaloReco_la_SOURCES = \
  RawClusterBuilderkV3.cc

libEICCaloReco_la_LIBADD = \
  -lCLHEP \
  -lphool \
  -lSubsysReco \
  -lgsl \
  -lgslcblas \
  -lg4vertex_io \
  -lcalo_io \
  -lphparameter

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libEICCaloReco.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
