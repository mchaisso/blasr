#include "algorithms/sorting/LightweightSuffixArray.h"
#include <iostream>
#include <string>
using namespace std;
#include "Types.h"
int main(int argc, char* argv[]) {
	UInt diffCover[] = {1,2,4};
	UInt diffCoverLength = 3;
	UInt diffCoverSize = 7;

	DiffCoverMu mu;
	UInt textLength = 26;
	mu.Initialize(diffCover, diffCoverLength, diffCoverSize, textLength);
	UInt dIndex = 0; // index in D-sample
	UInt tIndex = 0; // index in text
	UInt nDiffCover;
	nDiffCover = textLength / diffCoverSize + 1;
	UInt coverIndex;
	UInt d;
	bool done = false;
	
	for (coverIndex = 0; coverIndex < nDiffCover and done == false; coverIndex++) {
		for (d = 0; d < diffCoverLength and done == false; d++) {
			tIndex = coverIndex * diffCoverSize + diffCover[d];
			if (tIndex >= textLength) {
				done = true;
				cout << "done at " << coverIndex << " " << d << endl;
				break;
			}
			cout << "mu(" << dIndex << "): " << mu(tIndex) << endl;
			dIndex++;
		}
	}
}


								
