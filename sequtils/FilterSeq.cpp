#include <vector>
#include <string>
#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/NucConversion.h"
#include "../common/utils.h"
#include "../common/utils/StringUtils.h"


void PrintUsage() {
	cout << "usage: filterseq infile outfile [-toupper] [-line lineLength] [-minLength l]" << endl;
}



int main(int argc, char* argv[]) {

	string inFileName, outFileName;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	inFileName  = argv[1];
	outFileName = argv[2];
	int argi = 3;
	int doToUpper = 0;
	int lineLength = 50;
	DNALength minLength  = 0;
	vector<string> titlePatternsToRemove;
	while (argi < argc) {
		if (strcmp(argv[argi], "-toupper") == 0) {
			doToUpper = 1;
		}
		else if (strcmp(argv[argi], "-line") == 0) {
			lineLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minLength") == 0) {
			minLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-removeTitle") == 0) {
			titlePatternsToRemove.push_back(argv[++argi]);
		}
		else {
			PrintUsage();
			cout << "Bad option: " << argv[argi] << endl;
		}
		++argi;
	}

	FASTAReader reader;
	FASTASequence seq;
	reader.Init(inFileName);
	ofstream seqOut;
	CrucialOpen(outFileName, seqOut);
	int nBad = 0;
	while (reader.GetNext(seq)) {
		if (doToUpper) {
			seq.ToUpper();
		}
		if (seq.length < minLength) {
			continue;
		}
		int p;
		for (p = 0; p < titlePatternsToRemove.size(); p++) {
			if (ExactPatternMatch(seq.title, titlePatternsToRemove[p])) {
				continue;
			}
		}
		DNALength i;
		for (i = 0; i < seq.length; i++) { 
			seq.seq[i] = MaskedFourBit[seq.seq[i]];
			if (seq.seq[i] > 8) {
			++nBad;
			}
			seq.seq[i] = FourBitToAscii[seq.seq[i]];
		}
		seq.PrintSeq(seqOut, lineLength);
	}
}
