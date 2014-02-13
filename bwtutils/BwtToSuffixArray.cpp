#include "../common/datastructures/bwt/BWT.h"
#include "../common/datastructures/suffixarray/SuffixArray.h"
#include "../common/datastructures/suffixarray/SuffixArrayTypes.h"


#include <string>

using namespace std;
int main(int argc, char* argv[]) {
	
	string bwtFileName, saFileName;
	if (argc < 3) {
		cout << "usage: bwt2sa bwtfile safile " << endl;
		exit(1);
	}
	bwtFileName = argv[1];
	saFileName  = argv[2];

 	Bwt<PackedDNASequence, FASTASequence> bwt;	
	DNASuffixArray suffixArray;

	bwt.Read(bwtFileName);
	suffixArray.AllocateSuffixArray(bwt.bwtSequence.length-1);
	SAIndex index;
	for (index = 1; index < bwt.bwtSequence.length+1; index++) {
		suffixArray.index[index-1] = bwt.Locate(index);
	}
	suffixArray.Write(saFileName);
	
}
