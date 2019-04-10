#-----------------------------------------------------------------------------
#	Makefile for deploying the assembly cpp package to a given root location
#-----------------------------------------------------------------------------
SHELL = /bin/bash


INSTALL_LIB_DIR = $(PREFIX)/lib
INSTALL_BIN_DIR = $(PREFIX)/bin

# LIB_LIST = 

EXE_LIST =  \
  alignment \
  samutils

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

