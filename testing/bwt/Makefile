#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#

ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../../
endif

include ../../common.mk

all: bin make.dep testCreateFromSuffixArray testBuildOccBins testBuildPosTable testBWTCount testLFBacktrack testLocate

include ../../make.rules

include make.dep

testCreateFromSuffixArray: bin/testCreateFromSuffixArray
testBuildOccBins: bin/testBuildOccBins
testBuildPosTable: bin/testBuildPosTable
testBWTCount: bin/testBWTCount
testLFBacktrack: bin/testLFBacktrack
testLocate: bin/testLocate

bin/testCreateFromSuffixArray: bin/TestCreateFromSuffixArray.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testBuildOccBins: bin/TestBuildOccBins.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testBuildPosTable: bin/TestBuildPosTable.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testBWTCount: bin/TestBWTCount.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testLFBacktrack: bin/TestLFBacktrack.o
	$(CPP) $(CPPOPTS) $< -o $@


bin/testLocate: bin/TestLocate.o
	$(CPP) $(CPPOPTS) $< -o $@