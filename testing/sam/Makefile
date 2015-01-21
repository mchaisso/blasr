#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#

ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../../
endif

include ../../common.mk

all: bin make.dep testSAMReader

include ../../make.rules
include make.dep

PBCPP_DIR=../..


testSAMReader: bin/testSAMReader

bin/testSAMReader: bin/TestSAMReader.o
	$(CPP) $(CPPOPTS) $< -o $@
