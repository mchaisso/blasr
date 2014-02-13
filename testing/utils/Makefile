#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#
ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../../
endif

include ../../common.mk

all: bin make.dep testRandomSequence

include ../../make.rules
include make.dep


testRandomSequence: bin/testRandomSequence

bin/testRandomSequence: bin/TestRandomSequence.o
	$(CPP) $(CPPOPTS) $< -o $@

