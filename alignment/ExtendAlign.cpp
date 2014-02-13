#include "../common/algorithms/alignment/ExtendAlign.h"
#include "../common/datastructures/alignment/AlignmentCandidate.h"
#include "../common/files/ReaderAgglomerate.h"
#include "../common/FASTASequence.h"
#include "../common/algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "../common/algorithms/alignment/AlignmentPrinter.h"
#include <string>
#include <iostream>

int main(int argc, char* argv[]) {
	string fileAName, fileBName;
	if (argc < 3) {
		cout << "usage: extendAlign file1 fil2 [pos1 pos2] " << endl;
		exit(1);
	}

	fileAName = argv[1];
	fileBName = argv[2];
	int argi = 3;
	int aPos = 0;
	int bPos = 0;
	if (argc == 5) {
		aPos = atoi(argv[3]);
		bPos = atoi(argv[4]);
	}
	
	ReaderAgglomerate reader;
	reader.Initialize(fileAName);
	
	FASTASequence aSeq, bSeq;
	reader.GetNext(aSeq);
	reader.Initialize(fileBName);
	reader.GetNext(bSeq);
	
	DistanceMatrixScoreFunction<FASTASequence, FASTASequence> scoreFn;
	scoreFn.ins = 3;
	scoreFn.del = 3;
	scoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);

	vector<int>  scoreMat;
	vector<Arrow>pathMat;
	
	AlignmentCandidate<FASTASequence, FASTASequence> extendedAlignment;

	/*	ExtendAlignmentForward(aSeq, aPos,
												 bSeq, bPos,
												 5, //k
												 scoreMat, pathMat,
												 extendedAlignment,
												 scoreFn,
												 1, // don't bother attempting
												 // to extend the alignment
												 // if one of the sequences
												 // is less than 1 base long
												 2);

	extendedAlignment.qAlignedSeq.ReferenceSubstring(aSeq);
	extendedAlignment.tAlignedSeq.ReferenceSubstring(bSeq);

	//	extendedAlignment.qAlignedSeqPos = aPos;
	//	extendedAlignment.tAlignedSeqPos = bPos;

	StickPrintAlignment(extendedAlignment, aSeq, bSeq, cout);
	extendedAlignment.Clear();
	*/
	if (aPos == 0) { aPos = aSeq.length; }
	if (bPos == 0) { bPos = bSeq.length; }

	ExtendAlignmentReverse(aSeq, aPos,
												 bSeq, bPos,
												 5, //k
												 scoreMat, pathMat,
												 extendedAlignment,
												 scoreFn,
												 1, // don't bother attempting
												 // to extend the alignment
												 // if one of the sequences
												 // is less than 1 base long
												 2);

	extendedAlignment.qAlignedSeq.ReferenceSubstring(aSeq);
	extendedAlignment.tAlignedSeq.ReferenceSubstring(bSeq);

	//	extendedAlignment.qAlignedSeqPos = aPos;
	//	extendedAlignment.tAlignedSeqPos = bPos;

	StickPrintAlignment(extendedAlignment, aSeq, bSeq, cout);

	return 0;
}
