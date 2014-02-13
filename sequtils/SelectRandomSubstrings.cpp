#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/statistics/statutils.h"
#include "../common/utils.h"
#include <string>
#include <sstream>
using namespace std;
int main(int argc, char* argv[]) {
	if (argc < 5) {
		cout<< "usage: selectRandomSubstrings refFile outFile length num_seq" << endl;
		exit(1);
	}
	string refFileName = argv[1];
	string outFileName = argv[2];
	int    substrLength = atoi(argv[3]);
	int    numOutSeq   = atoi(argv[4]);
	vector<DNALength> cumulLengths;
	vector<FASTASequence*> refSequences;
	FASTAReader refReader;
	refReader.Initialize(refFileName);
	ofstream outFile;
	CrucialOpen(outFileName, outFile, std::ios::out);

	FASTASequence refSeq;
	cumulLengths.push_back(0);
	int seqIndex = 0;
	DNALength totalLength = 0;
	FASTASequence *seqPtr;
	while(refReader.GetNext(refSeq)) {
		seqPtr = new FASTASequence(refSeq);
		refSequences.push_back(seqPtr);
		totalLength = cumulLengths[seqIndex] + refSeq.length;
		cumulLengths.push_back(totalLength);
		seqIndex++;
	}
	int i;
	FASTASequence substr;
	if (substrLength > refSeq.length) {
		cout << "error,running with a substring length greater than the reference length." << endl;
		return 0;
	}
			
	for (i = 0; i < numOutSeq; i++) {

		int maxIts = 1000;
		bool printed = false;
		substr.length = substrLength;
		DNALength substrPos = RandomUnsignedInt(totalLength);
		while (!printed and maxIts > 0) {
			//
			// Determine which sequence this comes from.
			//
			for (seqIndex = 0; seqIndex < refSequences.size()-1; seqIndex++) {
				if (cumulLengths[seqIndex] <= substrPos and
						cumulLengths[seqIndex+1] > substrPos) {
					break;
				}
			}
			assert(seqIndex < refSequences.size());

			if (refSequences[seqIndex]->length < substrLength) {
				printed = false;
				continue;
			}
			assert(substrPos >= cumulLengths[seqIndex]);
			DNALength refSubstrPos = substrPos - cumulLengths[seqIndex];
			if (cumulLengths[seqIndex+1] - substrLength < substrPos) {
				refSubstrPos = cumulLengths[seqIndex+1] - substrLength;
			}
		  substr.seq = &refSequences[seqIndex]->seq[refSubstrPos];

		  int p;
		  for (p =0 ; p < substr.length; p++) {
		    if (substr[p] == 'N') {
					printed = false;
		      break;;
		    }
		  }
			
		  if (p == substr.length) {
				stringstream titlestrm;
				titlestrm << refSequences[seqIndex]->GetTitle() << "_" << refSubstrPos;
				substr.CopyTitle(titlestrm.str());
		    substr.PrintSeq(outFile);
		    printed = true;
		  }
		  --maxIts;
		}

	}

	return 0;
}
		

