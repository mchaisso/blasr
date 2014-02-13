#include "algorithms/sorting/LightweightSuffixArray.h"
#include <iostream>
#include <string>
using namespace std;
#include "Types.h"
int main(int argc, char* argv[]) {
	UInt D[] = {1,2,4};
	UInt lengthD = 3;
	UInt v = 7;

	UInt dh[7];

	BuildDiffCoverLookup(D, lengthD, v, dh);
	int i;
	for (i = 0; i < 7; i++) {
		cout << dh[i] << " "; 
	}
	cout << endl;
	return 0;
}
	
