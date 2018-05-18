
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#

ifeq ($(origin PREFIX), undefined)
  ifeq ($(origin SEYMOUR_HOME), undefined)
    PREFIX = /opt/seymour
  else
    PREFIX = $(SEYMOUR_HOME)
  endif
endif
ANALYSIS_HOME=$(PREFIX)/analysis

PBCPP_DIR = ..
SHELL = /bin/bash

#
# common.mk contains the configuration for this build setup.
#
include ../common.mk


#
# Define the targets before including the rules since the rules contains a target itself.
#

all: bin make.dep loadPulses  pls2fasta samtoh5 samtobas

include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif

#
# List of executables to deploy in the secondary analysis bunild
#



#
# build target for secondary build
# no point to building all if the secondary analysis
# system only deploys a few
# also keeps the build break frequency down
#
build-secondary: init $(DISTRIB_EXECS)
install-secondary: $(INSTALL_EXECS)


pls2fasta: bin/pls2fasta
loadPulses: bin/loadPulses
samtoh5: bin/samtoh5
samtobas: bin/samtobas




bin/pls2fasta: bin/PulseToFasta.o 
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/LoadPulses.o: LoadPulses.cpp 
	$(CPP) $(INCLUDEDIRS) $(CPPOPTS) -fno-trapping-math -g -c $<  -o $@

bin/loadPulses: bin/LoadPulses.o 
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/samtoh5: bin/SamToCmpH5.o
	$(CPP) -g $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/samtobas: bin/SamToBasH5.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)/

%-install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)/

clean:
	-rm -f bin/*
