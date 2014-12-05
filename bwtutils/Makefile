#
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

PBCPP_DIR = ../

#
# common.mk contains the configuration for this build setup.
#
include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

all: bin make.dep sa2bwt bwt2sa bwtLocateList bwtswPrint bwtpbPrint bwtquery

include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif

sa2bwt:      bin/sa2bwt
bwt2sa:      bin/bwt2sa
bwtLocateList: bin/bwtLocateList
bwtswPrint: bin/bwtswPrint
bwtpbPrint: bin/bwtpbPrint
bwtquery: bin/bwtquery

#
# List of executables to deploy in the secondary analysis build
#
INSTALL_EXES = sa2bwt bwt2sa
DISTRIB_EXES = $(addprefix bin/, $(INSTALL_EXES))
INSTALL_TARGETS := $(addsuffix -install, $(INSTALL_EXES))

install-secondary: $(INSTALL_TARGETS)

bin/sa2bwt: bin/SuffixArrayToBWT.o
	$(CPP) $(STATIC) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR)  $(HDF_REQ_LIBS)

bin/bwt2sa: bin/BwtToSuffixArray.o
	$(CPP) $(STATIC) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR)  $(HDF_REQ_LIBS)

bin/bwtLocateList: bin/BWTLocateList.o
	$(CPP) $(STATIC) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR)  $(HDF_REQ_LIBS)

bin/bwtswPrint: bin/PrintBWTSWBWT.o
	$(CPP) $(STATIC) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR)  $(HDF_REQ_LIBS)

bin/bwtpbPrint: bin/PrintPBBWT.o
	$(CPP) $(STATIC) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR)  $(HDF_REQ_LIBS)

bin/bwtquery: bin/BWTQuery.o
	$(CPP) $(STATIC) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR)  $(HDF_REQ_LIBS)

INSTALL_DIR ?= $(ANALYSIS_HOME)/bin

%-install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)/

install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 $(DISTRIB_EXES) $(INSTALL_DIR)/
