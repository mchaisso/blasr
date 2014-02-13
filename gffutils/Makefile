#
# common.mk contains the configuration for this build setup.
#

PBCPP_DIR = ../

include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

all: bin make.dep extractGFFSequences

include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif

extractGFFSequences:      bin/extractGFFSequences

bin/extractGFFSequences: bin/ExtractGFFSequences.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

