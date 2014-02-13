#include "algorithms/sorting/LightweightSuffixArray.h"
#include <iostream>
#include <string>
using namespace std;
#include "Types.h"
int main(int argc, char* argv[]) {
	UInt diffCover[] = {1,2,4};
	UInt diffCoverLength = 3;
	UInt diffCoverSize = 7;


	char text[] = "a rose is a rose is a rose";
	UInt textLength = strlen(text);
	UInt *index = new UInt[textLength+1];

	// The set d is given by 
	// Step 1.1 v-sort D-sample suffices
	UInt dIndex = 0; // index in D-sample
	UInt tIndex = 0; // index in text
	UInt nDiffCover;
	nDiffCover = textLength / diffCoverSize;
	UInt coverIndex;
	UInt d;
	bool done = false;
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
	UInt i;
	cout << "DSet indices: ";
	for (i = 0; i < dSetSize; i++) {
		cout << index[i] << " ";
	}
	UInt dSetSize = dIndex;
	MediankeyBoundedQuicksort(text, index, dIndex, 0, dSetSize, 0, diffCoverSize);

}
