# check that ROOTSYS is defined
ifndef ROOTSYS
  ROOTSYS = $(shell root-config --prefix)
  ifeq ($(ROOTSYS), )
    $(error running of root-config failed or reported null value)
  endif
endif

ROOTCFLAGS    = $(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTLIBS      = $(shell ${ROOTSYS}/bin/root-config --libs)
ROOTGLIBS     = $(shell ${ROOTSYS}/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -O2 -Wall -fPIC
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 
NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

.SUFFIXES: .cc,.C

# ================================================================================

all: fluxHistFits
#all: fitProg

fitProg: fitProg.o
	$(LD) $(LDFLAGS) -o fitProg fitProg.o $(GLIBS)

fluxHistFits: fluxHistFits.o
	$(LD) $(LDFLAGS) -o fluxHistFits fluxHistFits.o $(GLIBS)

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f fluxHistFits fitProg *.o
