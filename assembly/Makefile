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

PBCPP_DIR = ../

include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

EXECS = rm4ToRepeatGraph rm4ToPathList removeTransitiveOverlaps removeLowCoverageOverlaps rg2dot trimShortEnds printScaffolds untangle

# DISABLE for now
#cmpMatcher


all: bin make.dep $(EXECS)

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = rm4ToRepeatGraph
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))


include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif

rm4ToRepeatGraph:   bin/rm4ToRepeatGraph
rm4ToPathList:   bin/rm4ToPathList
removeTransitiveOverlaps: bin/removeTransitiveOverlaps
removeLowCoverageOverlaps: bin/removeLowCoverageOverlaps
trimShortEnds: bin/trimShortEnds
printScaffolds: bin/printScaffolds
untangle: bin/untangle

rg2dot: bin/rg2dot

bin/untangle: bin/Untangle.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/trimShortEnds: bin/TrimShortEnds.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/printScaffolds: bin/PrintScaffolds.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/rm4ToRepeatGraph: bin/RM4ToRepeatGraph.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/rm4ToPathList: bin/RM4ToPathList.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/removeTransitiveOverlaps: bin/RemoveTransitiveOverlaps.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/removeLowCoverageOverlaps: bin/RemoveLowCoverageOverlaps.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/rg2dot: bin/RGToDot.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

#
# Set up a default value for the install dir if one does
# not exist.
#
INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
BUILD_DIR ?= bin
install:
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)/

install-%:
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)/
