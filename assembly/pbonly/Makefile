#
# common.mk contains the configuration for this build setup.
#
#
ifeq ($(origin PREFIX), undefined)
  ifeq ($(origin SEYMOUR_HOME), undefined)
    PREFIX = /opt/seymour
  else
    PREFIX = $(SEYMOUR_HOME)
  endif
endif
ANALYSIS_HOME=$(PREFIX)/analysis

PBCPP_DIR = ../../

include ../../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

EXECS = buildEGraph

all: bin make.dep $(EXECS)

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = 
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))


include ../../make.rules

include make.dep

#testRGraph:   bin/testRGraph
buildRGraph: bin/buildRGraph
buildAGraph: bin/buildAGraph
buildEGraph: bin/buildEGraph


bin/buildEGraph: bin/BuildEBruijnGraph.o
	$(CPP) $(CPPOPTS) $< -o $@ 

bin/buildAGraph: bin/BuildABruijnGraph.o
	$(CPP) $(CPPOPTS) $< -o $@ 

bin/testRGraph: bin/TestRGraph.o
	$(CPP) $(CPPOPTS) $< -o $@ 

bin/buildRGraph: bin/BuildRGraph.o
	$(CPP) $(CPPOPTS) $< -o $@ 

#
# Set up a default value for the install dir if one does 
# not exist.
#
INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
BUILD_DIR ?= bin
install:
	/usr/bin/install -m 555 $(DISTRIB_EXECS) -t $(INSTALL_DIR) 

install-%:
	/usr/bin/install -m 555 bin/$* -t $(INSTALL_DIR)
