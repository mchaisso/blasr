#ifndef DATASTRUCTURES_BWT_POS_H_
#define DATASTRUCTURES_BWT_POS_H_

#include <fstream>
#include <vector>

#include "PackedHash.h"

#include "../../DNASequence.h"
#include "../../Types.h"
#include "../../utils/BitUtils.h"

using namespace std;
template< typename T_BWT_Sequence>
class Pos {
 public:
	static const unsigned int stride=8;
	PackedHash packedHash;
	vector<int> hashCount;
	vector<int> fullPos;
	int hasDebugInformation;
	void Write(ostream &out) {
		packedHash.Write(out);
	}

	void Read(istream &in) {
		packedHash.Read(in);
	}

	void InitializeFromSuffixArray(DNALength suffixArray[], DNALength suffixArrayLength) {
		DNALength p;
		packedHash.Allocate(suffixArrayLength);
		std::fill(hashCount.begin(), hashCount.end(), 0);
		for (p = 0; p < suffixArrayLength; p++ ){
			if (suffixArray[p] % stride==0){
				packedHash.AddValue(p,suffixArray[p]);
			}
		}
	}

	int Lookup(DNALength bwtPos, DNALength &seqPos) {
		return packedHash.LookupValue(bwtPos-1, seqPos);
	}

	

};


#endif
