#include <string>
#include <vector>
#include "../common/datastructures/matrix/FlatMatrix.h"
#include "../common/datastructures/alignment/Alignment.h"
#include "../common/algorithms/alignment/AlignmentUtils.h"
#include "../common/NucConversion.h"
#include "../common/utils.h"

int AlignmentStringToIndex(char *seq, int length) {
	int p;
	int index = 0;
	for (p = 0; p < length; p++ ) {
		if (seq[p] != '-') {
			index += TwoBit[seq[p]];
		}
		else {
			index +=4;
		}
		//
		// The equivalent of the shift in the tuple conversion case.
		//
		if (p <length-1)
			index *=5;
	}
	return index;
}


using namespace std;
int main(int argc, char* argv[]) {

	string alignmentFileName, tableOutName;
	int nNucs;
	if (argc < 4) {
		cout << "usage: tabulateAlignment compareSequenceAlignmentFile tableFile" << endl;
		exit(1);
	}
	alignmentFileName = argv[1];
	nNucs = atoi(argv[2]);
	tableOutName = argv[3];


	ifstream in;
	ofstream out;
	CrucialOpen(alignmentFileName, in);
	CrucialOpen(tableOutName, out, std::ios::out);

	int matDimSize = 1;
	int i;
	/*
	 Compute the length of one dimension of the matrix
	*/
	for (i = 0; i < nNucs; i++ )  {
		matDimSize = matDimSize * 5;
	}

	int a;
	int p;

	FlatMatrix2D<int> matrix(matDimSize + 1, matDimSize + 1);

	/*
	 * Tablulate the alignments one by one.
	 */
	
	CompSeqAlignment alignment;
	int alignLength;
	char *tSeq, *qSeq;
	while(ReadCompSeqAlignment(in, alignment)) {
		alignLength = alignment.tString.size();
		tSeq = (char*) alignment.tString.c_str();
		qSeq = (char*) alignment.qString.c_str();
		for (p = 0; p < alignLength - nNucs; p++ ) {
			//
			// Convert the 5-base nuc to an index.
			//
			int row, col;
			row = AlignmentStringToIndex(&tSeq[p], nNucs);
			col = AlignmentStringToIndex(&qSeq[p], nNucs);
			matrix.Set(row,col,matrix.Get(row,col) + 1);
		}
	}

	// add the column marginals to the matrix.
	int r, c;
	int rSum, cSum;
	for (r = 0; r < matDimSize; r++ ){ 
		cSum = 0;
		for (c = 0; c < matDimSize; c++ ){ 
			cSum += matrix.Get(r,c);
		}
		matrix.Set(r,matDimSize, cSum);
	}

	// add the row marginals to the matrix.
	
	for (c = 0; c < matDimSize; c++ ){ 
		rSum = 0;
		for (r = 0; r < matDimSize; r++ ){ 
			rSum += matrix.Get(r,c);
		}
		matrix.Set(matDimSize,c, rSum);
	}
	
	int tSum = 0;
	for (c = 0; c < matDimSize; c++ ){ 
		tSum += matrix.Get(matDimSize, c);
	}
	for (r = 0; r < matDimSize; r++) {
		tSum += matrix.Get(r,matDimSize);
	}
	matrix.Set(matDimSize, matDimSize, tSum);
	matrix.Print(out);
}

	
	
