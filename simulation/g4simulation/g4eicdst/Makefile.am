AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include

lib_LTLIBRARIES = \
   libg4eicdst.la

libg4eicdst_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lg4dst \
  -leiczdcbase \
  -leicpidbase

libg4eicdst_la_SOURCES = \
   g4eicdst.cc

g4eicdst.cc:
	echo "//*** this is a generated empty file. Do not commit, do not edit" > $@

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.cc \
  g4eicdst.cc

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libg4eicdst.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
