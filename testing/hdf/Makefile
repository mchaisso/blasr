#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#

ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../../
endif

include ../../common.mk

all: bin make.dep testHDFUtils testHDFAtom testHDFPlsReader testHDFBasReader testHDFRegionReader testHDFArrayWriter testHDF2DArrayWriter testHDFCmpWriter

include ../../make.rules
include make.dep


HDF5INCLUDEDIR = ../../../seymour/dist/common/include
HDF5LIBDIR = ../../../seymour/dist/common/lib/

testHDFUtils: bin/testHDFUtils
testHDFAtom: bin/testHDFAtom
testHDFPlsReader: bin/testHDFPlsReader
testHDFBasReader: bin/testHDFBasReader
#testHDFFile: bin/testHDFFile
testHDFRegionReader: bin/testHDFRegionReader
testHDFArrayWriter: bin/testHDFArrayWriter
testHDF2DArrayWriter: bin/testHDF2DArrayWriter
testHDFCmpWriter: bin/testHDFCmpWriter

#testCmpReader: bin/testCmpReader
#testHDFCCSReader: bin/testHDFCCSReader

bin/testHDFCmpWriter: bin/TestHDFCmpWriter.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIB) -l$(HDF5LIBCPP)

bin/testHDFReading: bin/TestHDFReading.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIB)

bin/testHDFArrayWriter: bin/TestHDFArrayWriter.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)

bin/testHDF2DArrayWriter: bin/TestHDF2DArrayWriter.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)

bin/cpptestHDFReading: bin/TestHDFReadingCPP.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)

bin/testHDFUtils: bin/TestHDFUtils.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)

bin/testHDFAtom: bin/TestHDFAtom.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)

bin/testHDFBasReader: bin/TestHDFBasReader.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)

bin/testHDFPlsReader: bin/TestHDFPlsReader.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)

bin/testHDFFile: bin/TestHDFFile.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB)

bin/testHDFRegionReader: bin/TestHDFRegionReader.o
	$(CPP) -static $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lpthread -lz

bin/testCmpReader: bin/TestHDFCmpReader.o
	$(CPP) -static $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lpthread -lz


bin/testHDFCCSReader: bin/TestHDFCCSReader.o
	$(CPP) -static $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lpthread -lz