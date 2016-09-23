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

PBCPP_DIR = ..

SHELL = /bin/bash		

#
# common.mk contains the configuration for this build setup.
#
include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

all: bin make.dep storeQualityByContext loadContextOut alchemy  sageify

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = storeQualityByContext loadContextOut alchemy
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))

include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif

storeQualityByContext:      bin/storeQualityByContext
loadContextOut: bin/loadContextOut
alchemy: bin/alchemy
sageify: bin/sageify

#
# List of executables to deploy in the secondary analysis build
# 

INSTALL_EXES = \
	storeQualityByContext \
	alchemy


DISTRIB_EXES = $(addprefix bin/, $(INSTALL_EXES))
#
# build target for secondary build
# no point to building all if the secondary analysis 
# system only deploys a few
# also keeps the build break frequency down
#
build-secondary: init $(DISTRIB_EXES)

INSTALL_TARGETS := $(addsuffix -install, $(INSTALL_EXES))

.PHONY: install-secondary

install-secondary: $(INSTALL_TARGETS)

bin/printLengthHistogram: bin/PrintLengthHistogramFromCmph5.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/storeQualityByContext: bin/StoreQualityByContext.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)  $(HDF_REQ_LIBS)

bin/alchemy: bin/Alchemy.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)  $(HDF_REQ_LIBS)

bin/loadContextOut: bin/LoadContextOut.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/sageify: bin/Sageify.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@


INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
BUILD_DIR ?= bin
install:
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)/

%-install:
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)/
install-%:
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)/
