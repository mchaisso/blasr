#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#
ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../../
endif

include ../../common.mk

all: bin make.dep testPackedDNA

include ../../make.rules
include make.dep

testPackedDNA: bin/testPackedDNA

bin/testPackedDNA: bin/TestPackedDNA.o
	$(CPP) $(CPPOPTS) $< -o $@