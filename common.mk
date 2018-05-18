#
# Definitions common to all make files.
#

HDF5INCLUDEDIR?=/home/mchaisso/software/include
HDF5LIBDIR?=/home/mchaisso/software/lib

INCLUDEDIRS = -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR) 

HDF5LIB    = hdf5
HDF5LIBCPP = hdf5_cpp
LINK_PROFILER = 
GCCOPTS = -O3 -Wno-div-by-zero $(INCLUDEDIRS) -fpermissive -mtune=native

HDF_REQ_LIBS= -lz -lpthread -ldl 
CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS) 
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  
CPP = g++ -std=c++11 -static
