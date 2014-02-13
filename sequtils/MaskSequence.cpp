#include "../common/files/AlignmentReader.h"
#include "../common/datastructures/alignment/AlignmentCandidate.h"
#include "../common/datastructures/metagenome/FASTATitleDictionary.h"
#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include <string>
typedef AlignmentCandidate<DNASequence,FASTQSequence> T_AlignmentCandidate;

int main(int argc, char* argv[]) {

	if (argc < 4) {
		cout << "usage maskSequence genomeFileName alignmentFileName maskedSeqName [-reverse] [-table]" << endl;
		cout << " -reverse Masks everything outside repeat regions." << endl;
		cout << " -table   Reads alignments as pairs of coordinates in a table rather than .out alignments." << endl;
		exit(1);
	}

	string genomeFileName, alignmentFileName, maskedSeqFileName;

	genomeFileName    = argv[1];
	alignmentFileName = argv[2];
	maskedSeqFileName = argv[3];
	bool reverse = false;
	bool table   = false;
	int argi = 4;
	while (argi < argc) {
		if (strcmp(argv[argi], "-reverse") == 0) {
			reverse = true;
		}
		else if (strcmp(argv[argi], "-table") == 0) {
			table = true;
		}
		else {
			cout << "ERROR, unknown option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}
	FASTAReader reader;
	AlignmentReader<T_AlignmentCandidate> alnReader;
	vector<FASTASequence> sequences;
	vector<T_AlignmentCandidate > alignments;
  FASTATitleDictionary titleDict;
	//
	// Initialize all files.
	//
	reader.Init(genomeFileName);

	ofstream maskedSequenceOut;
	CrucialOpen(maskedSeqFileName, maskedSequenceOut, std::ios::out);

	//
	// Grab data.
	//
	reader.ReadAllSequences(sequences);
	if (table == false) {
		alnReader.Init(alignmentFileName);
		alnReader.ReadAllAlignments(alignments);
	}
	else {
		ifstream tableIn;
		CrucialOpen(alignmentFileName, tableIn);
		while(tableIn) {
			T_AlignmentCandidate alignment;
			DNALength tEnd;
			if (!(tableIn >> alignment.tPos >> tEnd)) {
				break;
			}
			alignment.tLength = tEnd - alignment.tPos + 1;
			alignments.push_back(alignment);
		}
	}
	
	int s;
  for (s = 0; s < sequences.size(); s++) {
		titleDict.AddSequence(sequences[s]);
		sequences[s].ToUpper();
	}
	
	//
	// mask repeats.
	//
	int a;
	for (a = 0; a < alignments.size(); a++) {
		DNALength p;
		int seqIndex;
		seqIndex = -1;
		DNALength tStart, tEnd;
		tStart = tEnd = 0;

		if (table) {
			tStart = alignments[a].tPos;
			tEnd   = alignments[a].tPos + alignments[a].tLength;
			seqIndex = 0;
		}
		else if (titleDict.LookupSequence(alignments[a].tName, seqIndex)) {
			tStart = alignments[a].tPos;
			tEnd   = alignments[a].tPos + alignments[a].tLength;
		}
		else {
			cout << "Error! Aligned sequence " << seqIndex << " is not in the file." << endl;
		}
		if (reverse == false) {
			DNALength p;
			for (p = tStart; p < tEnd; p++ ) {
				sequences[seqIndex].seq[p] = 'N';
			}
		}
		else {
			cout << "setting to lower between " << tStart << " and " << tEnd << endl;
			for (p = tStart; p < tEnd; p++ ){
				sequences[0].seq[p] = tolower(sequences[0].seq[p]);
			}
		}
		
	}

	if (reverse) {
		DNALength p;
		for (p = 0; p < sequences[0].length; p++) {
			if (sequences[0].seq[p] >= 'A' and sequences[0].seq[p] <= 'Z') {
				sequences[0].seq[p] = 'N';
			}
			else {
				sequences[0].seq[p] = toupper(sequences[0].seq[p]);
			}
		}
	}


	for (s = 0; s < sequences.size(); s++ ) {
		sequences[s].PrintSeq(maskedSequenceOut);
	}

	

	return 0;
}


	

