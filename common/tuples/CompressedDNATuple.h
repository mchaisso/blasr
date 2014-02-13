#ifndef COMPRESSED_DNA_TUPLE_H_
#define COMPRESSED_DNA_TUPLE_H_

#include "DNATuple.h"
#include "../cmpseq/CompressedSequence.h"

template<typename T_Sequence>
class CompressedDNATuple : public DNATuple {
	static const unsigned char mask = 0xf;
 public:
	int FromStringLR(Nucleotide *strPtr, TupleMetrics &tm) {
		//
		// Make sure the sequence contains all valid characters.
		//

		if (!CompressedSequence<T_Sequence>::Only4BitACTG(strPtr, tm.tupleSize)) 
			return 0;
			
		if (tm.tupleSize == 0)
			return 1;

		tuple = 0;
		Nucleotide *p;
		Nucleotide *endPtr = &strPtr[tm.tupleSize - 1];
		for (p = strPtr; p != endPtr; p++) {
			tuple += TwoBit[*p & mask];
			tuple <<=2;
		}
		//
		// The tuple size is guaranteed to be at least 
		// 1, so it's safe to add the last value.
		// This cannot be in the previous loop since
		// the shift shouldn't happen.
		tuple += TwoBit[*p & mask];
		return 1;
	}


	int FromStringRL(Nucleotide *strPtr, TupleMetrics &tm) {
		
		if (!CompressedSequence<T_Sequence>::Only4BitACTG((CompressedNucleotide*)strPtr, tm.tupleSize))
			return 0;

		if (tm.tupleSize == 0)
			return 1;
			
		tuple = 0;
		Nucleotide *p;
		for (p = strPtr + tm.tupleSize - 1; p > strPtr; p--) {
			tuple += TwoBit[*p & mask];
			tuple <<=2;
		}
		//
		// The tuple size is guaranteed to be at least 
		// 1, so it's safe to add the last value.
		// This cannot be in the previous loop since
		// the shift shouldn't happen.
		tuple += TwoBit[*p & mask];
		return 1;
	}
};



#endif
