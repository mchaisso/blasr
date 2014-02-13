#include "algorithms/sorting/LightweightSuffixArray.h"
#include <iostream>
#include <string>
using namespace std;
#include "Types.h"
int main(int argc, char* argv[]) {
	UInt diffCover[] = {1,2,4};
	UInt diffCoverLength = 3;
	UInt diffCoverSize = 7;

	UInt textLength = 26;
	UInt dIndex = 0; // index in D-sample
	UInt tIndex = 0; // index in text
	UInt nDiffCover;
	nDiffCover = textLength / diffCoverSize + 1;
	UInt coverIndex;
	UInt d;
	bool done = false;
	
	char text[] = "a rose is a rose is a rose       ";
	//	UInt textLength = strlen(text);

	UInt *index = new UInt[textLength+1];

	for (coverIndex = 0; coverIndex < nDiffCover and done == false; coverIndex++) {
		for (d = 0; d < diffCoverLength and done == false; d++) {
			tIndex = coverIndex * diffCoverSize + diffCover[d];
			if (tIndex >= textLength) {
				done = true;
				break;
			}
			index[dIndex] = tIndex;
			dIndex++;
		}
	}
	UInt dSetSize = dIndex;
	MediankeyBoundedQuicksort((unsigned char*) text, index, dIndex, 0, dSetSize, 0, diffCoverSize);

	//
	// Step 1.2 Compute l^v(i) for all i \in D_n by traversing the
	// D-sample suffixes in lexicographic order and construct s^\prime
	// by setting s^\prime[\mu(i)] = l^v(i)
	//
	UInt *lexVNaming;
	lexVNaming = new UInt[dSetSize+1];
	DiffCoverMu mu;
	mu.Initialize(diffCover, diffCoverLength, diffCoverSize, textLength);

	UInt largestLexName;
	largestLexName = DiffCoverBuildLexNaming((unsigned char*) text, textLength,
																					 index, dSetSize, diffCover, diffCoverLength, diffCoverSize, 
																					 mu.diffCoverLookup,
																					 lexVNaming);

	cout << "lex " << diffCoverSize << " naming " << endl;
	UInt i;
	for (i = 0; i < dSetSize; i++) {
		cout << lexVNaming[i] << " ";
	}
	cout << endl;

}


								
