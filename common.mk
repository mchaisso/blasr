#
# Definitions common to all make files.
#

HDF5INCLUDEDIR ?= /usr/include
HDF5LIBDIR     ?= /usr/lib
HDF5INCLUDEDIR?=/usr/include
HDF5LIBDIR?=/usr/lib

INCLUDEDIRS = -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR) 

HDF5LIB    = hdf5
HDF5LIBCPP = hdf5_cpp
LINK_PROFILER = 
GCCOPTS = -std=c++11 -O3 -Wno-div-by-zero $(INCLUDEDIRS) -fpermissive -mtune=native

HDF_REQ_LIBS= -lz -lpthread -ldl 
CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS) 
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  
CPP = g++ -std=c++11 
