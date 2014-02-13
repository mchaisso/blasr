#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#
ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../../
endif

#
# common.mk contains the configuration for this build setup.
#
include ../../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

all: bin make.dep endAlignAllContigs pairAlign

include ../../make.rules

include make.dep

endAlignAllContigs:        bin/endAlignAllContigs
pairAlign: bin/pairAlign

bin/endAlignAllContigs: bin/EndAlignAllContigs.o
	$(CPP) $(CPPOPTS) $< -o $@  -L $(PBCPP_LIBDIR)

bin/pairAlign: bin/PairAlignAllShortContigs.o
	$(CPP) $(CPPOPTS) $< -o $@  -L $(PBCPP_LIBDIR)
  


clean:
	-rm bin/*
