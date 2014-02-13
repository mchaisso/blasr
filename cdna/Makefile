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

all: bin make.dep extractGenes gtfToGeneDB validateExons printDuplicates matchgtf

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = excrep evolve simpleShredder matchgtf
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))

include ../make.rules

include make.dep
extractGenes: bin/extractGenes
gtfToGeneDB: bin/gtfToGeneDB
validateExons: bin/validateExons
printDuplicates: bin/printDuplicates
matchgtf: bin/matchgtf

#printqual: bin/printqual

#
# List of executables to deploy in the secondary analysis build
# 

INSTALL_EXES = 

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

bin/extractGenes: bin/ExtractGenes.o
	$(CPP) $(CPPOPTS) $< -o $@  -lz -lpthread

bin/printDuplicates: bin/PrintDuplicates.o
	$(CPP) $(CPPOPTS) $< -o $@  -lz -lpthread

bin/validateExons: bin/ValidateExons.o
	$(CPP) $(CPPOPTS) $< -o $@  -lz -lpthread

bin/gtfToGeneDB: bin/GTFToGeneDB.o
	$(CPP) $(CPPOPTS) $< -o $@  -lz -lpthread

bin/matchgtf: bin/MatchGTFExons.o
	$(CPP) $(CPPOPTS) $< -o $@  -lz -lpthread

INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
BUILD_DIR ?= bin
install:
	/usr/bin/install -m 555 $(DISTRIB_EXECS) -t $(INSTALL_DIR) 

%-install:
	/usr/bin/install -m 555 bin/$* -t $(INSTALL_DIR)
install-%:
	/usr/bin/install -m 555 bin/$* -t $(INSTALL_DIR)
