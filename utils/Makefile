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

ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../
endif

include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

EXECS = join storeTuplePosList

all: bin make.dep $(EXECS)

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = join
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))

include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif

join:        bin/join
storeTuplePosList: bin/storeTuplePosList


bin/join: bin/Join.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -L $(PBCPP_LIBDIR)

bin/storeTuplePosList: bin/StoreTuplePosList.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -L $(PBCPP_LIBDIR)

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

clean:
	-rm -f bin/*
