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

PBCPP_DIR = ..

include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

EXECS =   blasrmc sawritermc blDotPlot

EXTENDED = guidedalign extendAlign

# DISABLE for now
#cmpMatcher


all: bin make.dep $(EXECS)

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = blasrmc sawritermc blDotPlot 
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))


include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif


blasrmc:        bin/blasrmc
blDotPlot:        bin/blDotPlot
sawritermc:           bin/sawritermc
STATIC=
bin/sawritermc: bin/SAWriter.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@


bin/blDotPlot: bin/DotPlot.o
	$(CPP) $(CPPOPTS) $^ $(STATIC) -o $@ 

bin/blasrmc: bin/Blasr.o
	$(CPP) $(STATIC) 	$(CPPOPTS) $< -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(LINK_PROFILER) $(HDF_REQ_LIBS) -lrt -o $@


#
# Set up a default value for the install dir if one does
# not exist.
#
INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
BUILD_DIR ?= bin
install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)

install-%:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)
