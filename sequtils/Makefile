#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#

ifeq ($(origin PREFIX), undefined)
  ifeq ($(origin SEYMOUR_HOME), undefined)
    PREFIX = /opt/seymour
  else
    PREFIX = $(SEYMOUR_HOME)
  endif
endif
ANALYSIS_HOME=$(PREFIX)/analysis

PBCPP_DIR = ..
SHELL = /bin/bash		

#
# common.mk contains the configuration for this build setup.
#
include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

all: bin make.dep catseq compressseq \
	filterseq bsdb catreads simpleShredder \
	countKMers evolve extractseq extractrep \
	printUngappedLengths dustoff splitContigs \
	printUnalignedSequences scramble maskSequence \
	selectRandomSubstrings printHPC \
	gcrank hammer checkfilt normgc excrep printTitleTable \
	makeBarcodeDatabase tileGenome filterBySubreadAgreement

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = excrep evolve simpleShredder
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))

include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif
printTitleTable: bin/printTitleTable
catreads:      bin/catreads
catseq:        bin/catseq
compressseq:      bin/compressseq
filterseq:     bin/filterseq
bsdb:           bin/bsdb
simpleShredder:bin/simpleShredder
countKMers: bin/countKmers
evolve: bin/evolve
extractseq: bin/extractseq
extractrep: bin/extractrep
#ourmatch: bin/ourmatch
printUngappedLengths: bin/printUngappedLengths
dustoff: bin/dustoff
splitContigs: bin/splitContigs
printUnalignedSequences: bin/printUnalignedSequences
scramble: bin/scramble
maskSequence: bin/maskSequence
selectRandomSubstrings: bin/selectRandomSubstrings
printHPC: bin/printHPC
gcrank: bin/gcrank
hammer: bin/hammer
checkfilt: bin/checkfilt
normgc: bin/normgc
excrep: bin/excrep
makeBarcodeDatabase: bin/makeBarcodeDatabase
tileGenome: bin/tileGenome
filterBySubreadAgreement: bin/filterBySubreadAgreement

#printqual: bin/printqual

#
# List of executables to deploy in the secondary analysis build
#

INSTALL_EXES = \
	bsdb \
	evolve

DISTRIB_EXES = $(addprefix bin/, $(INSTALL_EXES))
#
# build target for secondary build
# no point to building all if the secondary analysis
# system only deploys a few
# also keeps the build break frequency down
#
build-secondary: init $(DISTRIB_EXES)

INSTALL_TARGETS := $(addsuffix -install, $(INSTALL_EXES))

.PHONY: install-secondary

install-secondary: $(INSTALL_TARGETS)

bin/printTitleTable: bin/PrintTitleTable.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -lz -lpthread

bin/ourmatch: bin/OurMatch.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

bin/evolve: bin/Evolve.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/scramble: bin/ScrambleSequenceOrder.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/dustoff: bin/Dustoff.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/extractseq: bin/ExtractSeq.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/extractrep: bin/ExtractRepeats.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/printUngappedLengths: bin/PrintUngappedSequenceLengthHistogram.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/countKmers: bin/CountUniqueNMers.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/catseq: bin/ConcatenateSequences.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/catreads: bin/ConcatenateReads.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/compressseq: bin/Compress.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/filterseq: bin/FilterSeq.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/bsdb: bin/BuildSequenceDB.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/simpleShredder: bin/SimpleShredder.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/splitContigs: bin/SplitContigs.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/printUnalignedSequences: bin/PrintUnalignedSequences.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/maskSequence: bin/MaskSequence.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/selectRandomSubstrings: bin/SelectRandomSubstrings.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/printHPC: bin/PrintHomopolymerCount.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/gcrank: bin/PrintGCRank.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/hammer: bin/Hammer.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/checkfilt: bin/CheckForFilteredReads.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/normgc: bin/NormalizeGCContent.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/excrep: bin/ExciseRepeats.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/printrc: bin/ReverseComplement.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/makeBarcodeDatabase: bin/MakeBarcodeDatabase.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/tileGenome: bin/TileGenome.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/filterBySubreadAgreement: bin/FilterBySubreadAgreement.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -I $(HDF5INCLUDEDIR) -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(HDF_REQ_LIBS)

install-catseq:
	/usr/bin/install -m 555 bin/catseq $(INSTALL_DIR)/

install-compressseq:
	/usr/bin/install -m 555 bin/compresseq $(INSTALL_DIR)/

INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
BUILD_DIR ?= bin
install:
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)/

%-install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)
install-%:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 bin/$* $(INSTALL_DIR)
