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
		cout << "usage: sa2bwt genomeFileName suffixArray bwt [-debug]" << endl;
		exit(1);
	}
	string genomeFileName      = argv[1];
	string suffixArrayFileName = argv[2];
	string bwtFileName         = argv[3];
	int storeDebugInformation = 0;
	int argi = 4;
	while(argi < argc) {
		if (strcmp(argv[argi], "-debug") == 0) {
			storeDebugInformation = 1;
		}
		++argi;
	}
	
	ofstream bwtOutFile;
	CrucialOpen(bwtFileName, bwtOutFile, std::ios::out|std::ios::binary);
	
	FASTAReader reader;
	reader.Init(genomeFileName);
	FASTASequence seq;
	reader.ReadAllSequencesIntoOne(seq);
	
	
	DNASuffixArray suffixArray;
	suffixArray.Read(suffixArrayFileName);

	Bwt<PackedDNASequence, FASTASequence> bwt;
	bwt.InitializeFromSuffixArray(seq, suffixArray.index, storeDebugInformation ); 
	bwt.Write(bwtOutFile);

	return 0;
}

