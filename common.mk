#
# Definitions common to all make files.
#

HDF5INCLUDEDIR ?= /usr/include
HDF5LIBDIR     ?= /usr/lib
HDF5INCLUDEDIR=/net/eichler/vol5/home/mchaisso/software/hdf_memcheck/include/
HDF5LIBDIR=/net/eichler/vol5/home/mchaisso/software/hdf_memcheck/lib

INCLUDEDIRS = -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR) 

HDF5LIB    = hdf5
HDF5LIBCPP = hdf5_cpp
LINK_PROFILER = 
GCCOPTS =  -g -O3   -Wno-div-by-zero $(INCLUDEDIRS) -fpermissive

HDF_REQ_LIBS= -L$(HOME)/software/lib -lz -lpthread -ldl 
CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS) 
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  
CPP = g++
