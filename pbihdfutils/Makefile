
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

all: bin make.dep loadPulses writeHDFSubset pls2fasta extractReadsByTitle filterHDFPls pls2 testAdvance computeRefCoverage toAfg cmpView cmpInformation printmsa printvariants cmpH5Diff samtoh5 samtom4 samFilter samtobas writeTrimmedRegionTable

include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif

#
# List of executables to deploy in the secondary analysis bunild
#
DISTRIB_SET   = loadPulses writeHDFSubset pls2fasta toAfg samtoh5 samtom4 samFilter samtobas
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addsuffix -install, $(DISTRIB_SET))

#
# build target for secondary build
# no point to building all if the secondary analysis
# system only deploys a few
# also keeps the build break frequency down
#
build-secondary: init $(DISTRIB_EXECS)
install-secondary: $(INSTALL_EXECS)


writeHDFSubset: bin/writeHDFSubset
pls2: bin/pls2
pls2fasta: bin/pls2fasta
toAfg: bin/toAfg
filterHDFPls: bin/filterHDFPls
testAdvance: bin/testAdvance
loadPulses: bin/loadPulses
extractReadsByTitle: bin/extractReadsByTitle
pas: bin/pas
computeRefCoverage: bin/computeRefCoverage
cmpView: bin/cmpView
cmpInformation: bin/cmpInformation
printmsa: bin/printmsa
printvariants: bin/printvariants
cmpH5Diff: bin/cmpH5Diff
samtoh5: bin/samtoh5
samtom4: bin/samtom4
samtobas: bin/samtobas
samFilter: bin/samFilter
writeTrimmedRegionTable: bin/writeTrimmedRegionTable


bin/writeTrimmedRegionTable: bin/WriteTrimmedRegionTable.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/cmpView: bin/CmpH5View.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/cmpH5Diff: bin/DiffAlnInfoGroups.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/printmsa: bin/PrintMSA.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/printvariants: bin/PrintVariants.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/cmpInformation: bin/CmpH5Information.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/filterHDFPls: bin/FilterHDFPls.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/writeHDFSubset: bin/WriteHDFSubset.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testAdvance: bin/TestAdvance.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/pls2: bin/PrintHDFPlsSummary.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/toAfg: bin/ToAfg.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/pls2fasta: bin/PulseToFasta.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/LoadPulses.o: LoadPulses.cpp
	$(CPP) $(INCLUDEDIRS) $(CPPOPTS) -fno-trapping-math -g -c $<  -o $@

bin/loadPulses: bin/LoadPulses.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/pas: bin/PrintAlignmentStats.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/extractReadsByTitle: bin/ExtractReadsByTitle.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/computeRefCoverage: bin/ComputeRefCoverage.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/samtoh5: bin/SamToCmpH5.o
	$(CPP) -g $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/samtobas: bin/SamToBasH5.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)   -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/samtom4: bin/SamToM4.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)  $(HDF_REQ_LIBS)

bin/samFilter: bin/SamFilter.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@  -L$(HDF5LIBDIR)  $(HDF_REQ_LIBS)

INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)/

%-install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)/

clean:
	-rm -f bin/*
