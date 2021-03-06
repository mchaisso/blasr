#
# Definitions common to all make files.
#
# This should point to HDF v1.8.+

HDF5INCLUDEDIR?=/home/mchaisso/software/include
HDF5LIBDIR?=/home/mchaisso/software/lib

INCLUDEDIRS = -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR) 

HDF5LIB    = hdf5
HDF5LIBCPP = hdf5_cpp
LINK_PROFILER = 
GCCOPTS ?= -O3 -std=c++11

HDF_REQ_LIBS= -lz -lpthread -ldl 
CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS) 
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  
CPP = g++ 
