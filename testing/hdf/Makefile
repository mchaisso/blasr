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
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIB) -l$(HDF5LIBCPP) $(HDF_REQ_LIBS)

bin/testHDFReading: bin/TestHDFReading.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testHDFArrayWriter: bin/TestHDFArrayWriter.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testHDF2DArrayWriter: bin/TestHDF2DArrayWriter.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/cpptestHDFReading: bin/TestHDFReadingCPP.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testHDFUtils: bin/TestHDFUtils.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testHDFAtom: bin/TestHDFAtom.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testHDFBasReader: bin/TestHDFBasReader.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testHDFPlsReader: bin/TestHDFPlsReader.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testHDFFile: bin/TestHDFFile.o
	$(CPP) $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testHDFRegionReader: bin/TestHDFRegionReader.o
	$(CPP) -static $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/testCmpReader: bin/TestHDFCmpReader.o
	$(CPP) -static $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)


bin/testHDFCCSReader: bin/TestHDFCCSReader.o
	$(CPP) -static $(CPPOPTS) $< -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)
