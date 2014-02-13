#include "FASTASequence.h"
#include "FASTAReader.h"

#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "datastructures/suffixarray/SuffixArray.h"
#include "datastructures/bwt/BWT.h"
#include "utils.h"

#include <iostream>
#include <fstream>


using namespace std;
int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: bwtLocateList bwtName querySeqFile" << endl;
		exit(1);
	}
	string bwtFileName      = argv[1];
	string querySeqFileName = argv[2];
	bool doPrintResults = false;
	int maxCount = 0;
	int argi = 3;
	bool countOnly = false;
	while(argi < argc) {
		if (strcmp(argv[argi], "-print") == 0) {
			doPrintResults = true;
		}
		else if (strcmp(argv[argi], "-max") == 0) {
			maxCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-count") == 0) {
			countOnly = true;
		}
		else {
			cout << "bad option: " << argv[argi] << endl;
		}
		++argi;
	}

 	Bwt<PackedDNASequence, FASTASequence> bwt;
	bwt.Read(bwtFileName);

	FASTAReader queryReader;
	queryReader.Init(querySeqFileName);
	FASTASequence seq;
	int seqIndex = 0;
	vector<DNALength> positions;
	while(queryReader.GetNext(seq)) {
		positions.clear();
		if (countOnly == false) {
			bwt.Locate(seq, positions, maxCount);
		}
		else {
			DNALength sp,ep;
			bwt.Count(seq, sp, ep);
		}
		//		cout << "matched " << positions.size() << " positions." << endl;
		if (doPrintResults) {
			int i;
			for (i = 0; i < positions.size(); i++ ){
				cout << positions[i] << " ";
			}
			cout << endl;
		}
		++seqIndex;
	}
		//	float wordCountsPerLookup = (bwt.bwtSequence.nCountInWord *1.0) / bwt.bwtSequence.nCountNuc;
		//	cout << "word counts per lookup: " << wordCountsPerLookup << endl;
	return 0;
}


