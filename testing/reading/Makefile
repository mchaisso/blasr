#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#
ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../../
endif

include ../../common.mk

all: bin make.dep testra testCmpSeqAlignmentReader  testFASTAReader testRandomSequence

include ../../make.rules
include make.dep

testra: bin/testra bin/testCmpSeqAlignmentReader
testCmpSeqAlignmentReader: bin/testCmpSeqAlignmentReader
readTitleTable: bin/readTitleTable
testPSLReader: bin/testPSLReader
testFASTAReader: bin/testFASTAReader


bin/testFASTAReader: bin/TestFASTAReader.o
	$(CPP) $(CPPOPTS) $< -o $@ -fPIC

bin/readTitleTable: bin/ReadTitleTable.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testra: bin/TestReaderAgglomerate.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIB) -l$(HDF5LIBCPP) $(HDF_REQ_LIBS)

bin/testCmpSeqAlignmentReader: bin/TestCmpSeqAlignmentReader.o
	$(CPP) $(CPPOPTS) $< -o $@

bin/testPSLReader: bin/TestPSLReader.o
	$(CPP) $(CPPOPTS) $< -o $@
