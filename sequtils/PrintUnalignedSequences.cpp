#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/utils.h"

#include <string>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "usage: printUnalignedSequences in.fa alignfile outfile" <<endl;
		exit(1);
	}
	string inFileName, outFileName, alignmentFileName;
	inFileName = argv[1];
	alignmentFileName = argv[2];
	outFileName = argv[3];

	FASTAReader reader;
	reader.Init(inFileName);
	FASTASequence genome;

	reader.GetNext(genome);

	vector<bool> alignedSeq;
	alignedSeq.resize(genome.length);
	std::fill(alignedSeq.begin(), alignedSeq.end(), false);
	ifstream alignmentIn;
	ofstream gapFileOut;
	CrucialOpen(alignmentFileName, alignmentIn, std::ios::in);
	CrucialOpen(outFileName, gapFileOut, std::ios::out);
	string line;
	string s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12;
	DNALength genomeStart, genomeEnd;

	// mask off all portions of the genome that are not contained in an alignment
	while(alignmentIn) {
		alignmentIn >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9 >> genomeStart >> genomeEnd >> s12;
		DNALength readPos;
		for (readPos = genomeStart; readPos <= genomeEnd; readPos++ ){
			alignedSeq[readPos] = true;
		}
	}

	FASTASequence gapSeq;
	DNALength gapStart, gapEnd;
	stringstream titleStrm;
	int gapIndex = 0;
	int gapLength = 0;
	
	DNALength genomePos = 0;
	while(genomePos < genome.length) { 
		gapEnd = gapStart = genomePos;
		while (gapStart < genome.length and alignedSeq[gapStart] == true) { gapStart++;}
		gapEnd = gapStart + 1;
		while (gapEnd < genome.length and alignedSeq[gapEnd] == false) { gapEnd++;}
		titleStrm.str("");
		titleStrm << gapIndex << "_" << gapEnd - gapStart << "_" << gapStart << "_" << gapEnd;
		++gapIndex;
		if (gapStart >= genome.length) 
			break;
		if (gapEnd - gapStart > 0) {
			gapSeq.CopyTitle(titleStrm.str());
			gapSeq.seq = &genome.seq[gapStart];
			gapSeq.length = gapEnd - gapStart;
			gapSeq.PrintSeq(gapFileOut);
			genomePos = gapEnd;
			cerr << gapEnd << " " << gapStart << " " << gapEnd - gapStart << " " << genomePos << endl;
		}
		else {
			genomePos++;
		}
	}
}
			

			
	
	
