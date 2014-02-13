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

EXECS = testSuffixArray \
	testRadixSort \
	testSequencedb \
	testMultikeyQuicksort \
	testBuildDiffCoverLookup \
	testDiffCoverDelta \
	testlsa \
	testdsort \
	testmu

all: bin make.dep $(EXECS)

BUILTEXECS = $(addprefix bin/, $(EXECS))

include ../../make.rules

include make.dep

testSuffixArray:    bin/testSuffixArray
testRadixSort:      bin/testRadixSort
testSequencedb:     bin/testSequencedb
testMultikeyQuicksort: bin/testMultikeyQuicksort
testBuildDiffCoverLookup: bin/testBuildDiffCoverLookup
testDiffCoverDelta: bin/testDiffCoverDelta
testdsort: bin/testdsort
testlsa: bin/testlsa
testmu: bin/testmu
testlexnaming: bin/testlexnaming
testlexorder: bin/testlexorder

CPPOPTS = -g
bin/testSuffixArray: bin/TestSuffixArray.o $(PBCPP_DIR)/common/algorithms/sorting/qsufsort.o
	$(CPP) $(CPPOPTS) $< -o $@ $(PBCPP_DIR)/common/algorithms/sorting/qsufsort.o

bin/testMultikeyQuicksort: bin/TestMultikeyQuicksort.o 
	$(CPP) $(CPPOPTS) $< -o $@

bin/testRadixSort: bin/TestRadixSort.o
	$(CPP) $(CPPOPTS) $< -o $@ 

bin/testSequencedb: bin/TestSequenceDB.o
	$(CPP) $(CPPOPTS) $< -o $@ 

bin/testBuildDiffCoverLookup: bin/TestBuildDiffCoverLookup.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testDiffCoverDelta: bin/TestDiffCoverDelta.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testlsa: bin/TestLightweightSuffixArray.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testdsort: bin/TestDSetSort.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testmu: bin/TestDiffCoverMu.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testlexnaming: bin/TestLexNaming.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testlexorder: bin/TestLexOrder.o
	$(CPP) $(CPPOPTS) $< -o $@
