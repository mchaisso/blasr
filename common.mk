#
# Definitions common to all make files.
#

<<<<<<< HEAD
HDF5INCLUDEDIR ?= /usr/include
HDF5LIBDIR     ?= /usr/lib
HDF5INCLUDEDIR?=/net/eichler/vol5/home/mchaisso/software/hdf/include/
HDF5LIBDIR?=/net/eichler/vol5/home/mchaisso/software/hdf/lib
=======
HDF5INCLUDEDIR?=/home/mchaisso/software/include
HDF5LIBDIR?=/home/mchaisso/software/lib
>>>>>>> 3572cfb26594ed07750e1781e6ef9f2ff4356583

INCLUDEDIRS = -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR) 

HDF5LIB    = hdf5
HDF5LIBCPP = hdf5_cpp
LINK_PROFILER = 
GCCOPTS = -std=c++11 -O3 -Wno-div-by-zero $(INCLUDEDIRS) -fpermissive -mtune=native

HDF_REQ_LIBS= -lz -lpthread -ldl 
CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS) 
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  
CPP = g++ -std=c++11 
