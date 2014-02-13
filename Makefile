#-----------------------------------------------------------------------------
#	Makefile for deploying the assembly cpp package to a given root location
#-----------------------------------------------------------------------------
SHELL = /bin/bash

ifeq ($(origin ASSEMBLY_HOME), undefined)
PREFIX = $(SEYMOUR_HOME)/analysis
else
PREFIX = $(ASSEMBLY_HOME)
endif

INSTALL_LIB_DIR = $(PREFIX)/lib
INSTALL_BIN_DIR = $(PREFIX)/bin

# LIB_LIST = 

EXE_LIST =  \
  alignment \
  sequtils \
	pbihdfutils \
	bwtutils \
  simulator \
  samutils

# Not currently used - AAK
# 	pelusa \
#	jabon \
# All common files should be headers- AAK
#  common \


BUILT_EXES := $(foreach subdir, $(EXE_LIST), $(subdir)/build/$(subdir))

#--- Building
BUILD_TARGETS := $(addsuffix -build, $(EXE_LIST))
%-build: 
	@[[ -e $*/Makefile ]] && $(MAKE) -C $* -f Makefile UP=../ || true
build: $(BUILD_TARGETS) 

#--- Installing --- TODO put in something for shared libraries
INSTALL_TARGETS := $(addsuffix -install, $(EXE_LIST))
%-install:
	@[[ -e $*/Makefile ]] && $(MAKE) -C $* -f Makefile install INSTALL_DIR=$(INSTALL_BIN_DIR) || true
install: $(INSTALL_TARGETS)

cramtests:
	cram --shell=/bin/bash ctest/*.t

#--- Cleaning
CLEAN_TARGETS := $(addsuffix -clean, $(EXE_LIST))
%-clean:
	@[[ -e $*/Makefile ]] && $(MAKE) -C $* -f Makefile clean || true
clean: $(CLEAN_TARGETS)

