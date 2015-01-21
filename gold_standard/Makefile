#
# common.mk contains the configuration for this build setup.
#

ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../
endif

include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

EXECS = exhalign

all: bin make.dep $(EXECS)

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = exhalign
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))

include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif



exhalign:        bin/exhalign

bin/exhalign: bin/ExhaustiveAlign.o
	$(CPP) -O3 $< $(STATIC) -o $@  -L $(HDF5LIBDIR) -l $(HDF5LIBCPP)  -l $(HDF5LIB) $(HDF_REQ_LIBS) 

#
# Set up a default value for the install dir if one does 
# not exist.
#
INSTALL_DIR ?= $(ASSEMBLY_HOME)/bin
BUILD_DIR ?= bin
install:
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)/

install-%:
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)/

clean:
	-rm -f bin/*
