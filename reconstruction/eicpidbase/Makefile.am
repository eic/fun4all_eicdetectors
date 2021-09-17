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
	EICPIDDefs.h \
	EICPIDParticle.h \
	EICPIDParticleContainer.h \
	EICPIDParticlev1.h 

lib_LTLIBRARIES = \
  libeicpidbase.la

libeicpidbase_la_SOURCES = \ \
  $(ROOTDICTS) \
	EICPIDDefs.cc \
	EICPIDParticle.cc \
	EICPIDParticleContainer.cc \
	EICPIDParticlev1.cc

libeicpidbase_la_LIBADD = \
  -lphool  


ROOTDICTS = \
  EICPIDParticle_Dict.cc \
  EICPIDParticleContainer_Dict.cc \
  EICPIDParticlev1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  EICPIDParticle_Dict_rdict.pcm \
  EICPIDParticleContainer_Dict_rdict.pcm \
  EICPIDParticlev1_Dict_rdict.pcm


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libeicpidbase.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
