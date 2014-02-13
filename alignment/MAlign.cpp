#include <assert.h>
#include <string>
#include <iostream>
#include <vector>

#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/algorithms/alignment.h"
#include "../common/algorithms/alignment/BellEndAlign.h"
#include "../common/datastructures/matrix/FlatMatrix.h"


int main(int argc, char* argv[]) {
	
	string xyzIn;
	if (argc <= 1) {
		cout << "usage: malign reads.fasta.  This file should have 3 sequences.  The first is " <<endl
				 << "       a prefix of the second, and the third is a suffix of the second." << endl << endl;
		exit(1);
	}
	xyzIn = argv[1];
	FASTAReader reader(xyzIn);

	FASTASequence seqX, seqY, seqZ;

	reader.GetNext(seqX);
	reader.GetNext(seqY);
	reader.GetNext(seqZ);
	
	seqX.ToThreeBit();
	seqY.ToThreeBit();
	seqZ.ToThreeBit();
	//	Alignment alignment;
	vector<Arrow> optAlignment;
	FlatMatrix3D<int> scoreMat;
	FlatMatrix3D<Arrow> pathMat;

	BellEndAlign(seqX, seqY, seqZ, SMRTDistanceMatrix, 10, scoreMat, pathMat, optAlignment);

	string stringX, stringY, stringZ;

	CreateBellEndAlignmentStrings(seqX, seqY, seqZ, optAlignment, 
															 stringX, stringY, stringZ);
	
	PrintMAlignStrings(stringX, stringY, stringZ, cout);
}
