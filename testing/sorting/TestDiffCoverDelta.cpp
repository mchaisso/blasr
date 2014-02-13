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

	DiffCoverDelta delta;
	delta.Initialize(D, lengthD, v);

	cout << "delta(2,12) " << delta(2,12) << " should be 6" << endl;
	return 0;
}
	
