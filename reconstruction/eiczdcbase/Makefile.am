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
	RawTowerZDC.h \
	RawTowerZDCv1.h \
	RawTowerZDCDefs.h \
	RawTowerZDCContainer.h \
	RawTowerZDCGeom.h \
	RawTowerZDCGeomv1.h \
	RawTowerZDCGeomContainer.h \
	RawTowerZDCDeadMap.h

lib_LTLIBRARIES = \
  libeiczdcbase.la

libeiczdcbase_la_SOURCES = \ \
  $(ROOTDICTS) \
	RawTowerZDC.cc \
	RawTowerZDCv1.cc \
	RawTowerZDCContainer.cc \
	RawTowerZDCGeom.cc \
	RawTowerZDCGeomv1.cc \
	RawTowerZDCGeomContainer.cc \
	RawTowerZDCDeadMap.cc

libeiczdcbase_la_LIBADD = \
  -lphool  

ROOTDICTS = \
	RawTowerZDC_Dict.cc \
	RawTowerZDCv1_Dict.cc \
	RawTowerZDCContainer_Dict.cc \
	RawTowerZDCGeom_Dict.cc \
	RawTowerZDCGeomv1_Dict.cc \
	RawTowerZDCGeomContainer_Dict.cc \
	RawTowerZDCDeadMap_Dict.cc
pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
	RawTowerZDC_Dict_rdict.pcm \
	RawTowerZDCv1_Dict_rdict.pcm \
	RawTowerZDCContainer_Dict_rdict.pcm \
	RawTowerZDCGeom_Dict_rdict.pcm \
	RawTowerZDCGeomv1_Dict_rdict.pcm \
	RawTowerZDCGeomContainer_Dict_rdict.pcm  \
	RawTowerZDCDeadMap_Dict_rdict.pcm

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libeiczdcbase.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
