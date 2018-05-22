
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

all: bin make.dep printGaps samtom svrealigner

include ../make.rules

include make.dep

#
# List of executables to deploy in the secondary analysis build
# 
DISTRIB_SET   = printGaps
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

printGaps: bin/printGaps

printGappyRegions: bin/printGappyRegions
printMaxAlignedRead: bin/printMaxAlignedRead
printCoverage: bin/printCoverage
samtom: bin/samtom
svrealigner: bin/svrealigner

bin/printGaps: bin/PrintGaps.o
	$(CPP) $(CPPOPTS) $< -o $@ -lpthread 

bin/samtom: bin/SamToM0.o
	$(CPP) $(CPPOPTS) $< -o $@ -lpthread

bin/svrealigner: bin/SVRealigner.o
	$(CPP) $(CPPOPTS) $< -o $@ -lpthread

INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
install:
	/usr/bin/install -m 555 $(DISTRIB_EXECS) -t $(INSTALL_DIR) 

%-install:
	/usr/bin/install -m 555 bin/$* -t $(INSTALL_DIR)

clean:
	-rm bin/*
