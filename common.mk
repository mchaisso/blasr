#
# Definitions common to all make files.
#

HDF5INCLUDEDIR ?= /usr/include
HDF5LIBDIR     ?= /usr/lib

INCLUDEDIRS = -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR) 

HDF5LIB    = hdf5
HDF5LIBCPP = hdf5_cpp
LINK_PROFILER = 
GCCOPTS = -O3 -Wno-div-by-zero $(INCLUDEDIRS) -fpermissive  -I$(HOME)/software/include 
#
HDF_REQ_LIBS= -L$(HOME)/software/lib -lz -lpthread -ldl 
CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS) 
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  
CPP = g++
