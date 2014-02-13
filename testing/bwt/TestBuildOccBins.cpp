#include "FASTASequence.h"
#include "FASTAReader.h"

#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "datastructures/suffixarray/SuffixArray.h"
#include "datastructures/bwt/BWT.h"

#include <iostream>
#include <fstream>


using namespace std;
int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: testBuildOccBins genomeFileName suffixArray" << endl;
		exit(0);
	}
	string genomeFileName      = argv[1];
	string suffixArrayFileName = argv[2];
	
	FASTAReader reader;
	reader.Init(genomeFileName);
	FASTASequence seq;
	reader.GetNext(seq);
	

	DNASuffixArray suffixArray;
	
	suffixArray.Read(suffixArrayFileName);

	Bwt<PackedDNASequence, FASTASequence> bwt;
	//bwt.InitializeFromSuffixArray(seq, suffixArray.index);

	bwt.InitializeBWTStringFromSuffixArray(seq, suffixArray.index);
	bwt.occ.Initialize(bwt.bwtSequence, 4096, 64);
	bwt.occ.PrintBins(cout);
}


