#include "algorithms/sorting/MultikeyQuicksort.h"
#include "FASTAReader.h"
#include "FASTASequence.h"

#include <vector>
#include <string>
using namespace std;


int main(int argc, char* argv[]) {
	string seqFileName;
	int    bound;
	if (argc < 3) {
		cout << "usage: testMultikeyQuicksort seq bound" << endl;
		exit(0);
	}

	seqFileName = argv[1];
	bound = atoi(argv[2]);

	FASTAReader reader;
	reader.Init(seqFileName);
	reader.SetSpacePadding(bound);
	FASTASequence seq;

	reader.GetNext(seq);
	TransformSequenceForSorting(seq.seq, seq.length, bound);

	UInt *index = new UInt[seq.length];
	UInt i;
	for (i = 0; i < seq.length; i++) {
		index[i] = i;
	}
	MediankeyBoundedQuicksort(seq.seq, index, seq.length,0, seq.length, 0, bound);

	TransformBackSequence(seq.seq, seq.length);
	/*
	// print all suffices of the text
	UInt sufIndex;
	for (sufIndex = 0; sufIndex < seq.length; sufIndex++ ){
		DNASequence sufSeq;
		sufSeq.seq = &seq.seq[index[sufIndex]];
		sufSeq.length = seq.length - index[sufIndex];
		sufSeq.PrintSeq(cout);
		}*/

}

	
	
	
