#include "FASTASequence.h"
#include "FASTAReader.h"

#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "datastructures/suffixarray/SuffixArray.h"
#include "datastructures/bwt/BWT.h"

#include <iostream>
#include <fstream>


using namespace std;
int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "usage: testBWTCount genomeFileName suffixArray querySequence" << endl;
		exit(0);
	}
	string genomeFileName      = argv[1];
	string suffixArrayFileName = argv[2];
	string queryFileName       = argv[3];

	FASTAReader reader;
	reader.Init(genomeFileName);
	FASTASequence seq;
	reader.GetNext(seq);
	
	FASTAReader queryReader;
	queryReader.Initialize(queryFileName);
	FASTASequence query;
	queryReader.GetNext(query);
	
	DNASuffixArray suffixArray;
	
	suffixArray.Read(suffixArrayFileName);

	Bwt<PackedDNASequence, FASTASequence> bwt;
	//
	// Build the whole thing.
	//
	bwt.InitializeFromSuffixArray(seq, suffixArray.index);
	bwt.bwtSequence.PrintUnpacked(cout);
	cout <<endl;
	int nocc = bwt.Count(query);
	cout << nocc << endl;
}


