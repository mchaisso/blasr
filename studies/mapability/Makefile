#
# Configure the base directory fo the secondary c++ source, if it is
# not already specified.
#
ifeq ($(origin PBCPP_DIR), undefined)
PBCPP_DIR = ../../
endif

#
# common.mk contains the configuration for this build setup.
#
include ../../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#
CCOPTS = -O3 
all: bin make.dep countKmers nklMapability sequenceToTupleList countWordUniqueness printWordCount countKmerMultiplicity

include ../../make.rules

include make.dep

countKmers:        bin/countKmers
#extractRepeat: bin/extractRepeat
nklMapability: bin/nklMapability
sequenceToTupleList: bin/sequenceToTupleList
countWordUniqueness: bin/countWordUniqueness
printWordCount: bin/printWordCount

bin/countKmers: bin/CountUniqueKMersWithSuffixArray.o
	$(CPP) $(CPPOPTS) $< -o $@  -L $(PBCPP_LIBDIR)

bin/countWordUniqueness: bin/CountWordUniqueness.o
	$(CPP) $(CPPOPTS) $< -o $@  -L $(PBCPP_LIBDIR) -lpthread

bin/sequenceToTupleList: bin/SequenceToTupleList.o
	$(CPP) $(CPPOPTS) $< -o $@  -L $(PBCPP_LIBDIR)

bin/extractRepeat: bin/ExtractRepeat.o
	$(CPP) $(CPPOPTS) $< -o $@  

bin/nklMapability: bin/ComputeNKLMapability.o
	$(CPP) $(CPPOPTS) $< -o $@ -lpthread

bin/printWordCount: bin/PrintWordCount.o
	$(CPP) $(CPPOPTS) $< -o $@ -lpthread

clean:
	-rm bin/*