#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = dEdxConversionDict
SHLIB = libdEdxConversion.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

CXX=G++
CXXFLAGS= -g
# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += $(shell root-config --cflags)
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += -I$(LARLITE_COREDIR)
INCFLAGS += -I$(LARLITE_USERDEVDIR)
INCFLAGS += -I$(LARCV_APPDIR)/UBWireTool

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
LDFLAGS += -L.
LDFLAGS += $(shell root-config --libs)
LDFLAGS += -L$(LARLITE_LIBDIR)
LDFLAGS += $(shell larlite-config --libs)
include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE


test: test.cxx
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $^ -o test Cal_Test.o $(LDFLAGS)

clntst: 
	rm -rf test.dSYM
	rm test


