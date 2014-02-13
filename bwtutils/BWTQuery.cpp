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
		cout << "usage: bwtQuery bwtName querySeqFile" << endl;
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
	vector<DNALength> spv, epv;

	while(queryReader.GetNext(seq)) {
		positions.clear();

    epv.clear();
    spv.clear();
    bwt.Count(seq, spv, epv);
    if (spv.size() > 1 and epv.size() > 1) {
      assert(spv.size() == epv.size());
      int matchLength = spv.size();
      int last = matchLength - 1;
      if (epv[last] - spv[last] < maxCount) {
        bwt.Locate(spv[last], epv[last], positions);
        cout <<  matchLength << " " << epv[last] - spv[last] + 1 << endl;
      }
    }
	}
	return 0;
}


