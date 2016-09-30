#ifndef DNA_TUPLE_H_
#define DNA_TUPLE_H_


#include <vector>

#include "TupleList.h"
#include "BaseTuple.h"
#include "TupleOperations.h"

#include "Types.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "NucConversion.h"
using namespace std;


class DNATuple : public BaseTuple {
 public:
	DNALength pos;
 DNATuple() : BaseTuple() {
		pos = 0;
	}
	int FromStringLR(Nucleotide *strPtr, TupleMetrics &tm) {
		DNASequence tmpSeq;
		tmpSeq.seq = strPtr;
		tmpSeq.length = tm.tupleSize;
		if (!OnlyACTG(tmpSeq))
			return 0;

		if (tm.tupleSize == 0)
			return 1;

		tuple = 0;
		Nucleotide *p;
		Nucleotide *endPtr = &strPtr[tm.tupleSize - 1];
		for (p = strPtr; p != endPtr; p++) {
			// If it is not possible to convert this string, return null.
			if (ThreeBit[*p] > 3) {
				return 0;
			}
			tuple += TwoBit[*p];
			tuple <<=2;
		}
		//
		// The tuple size is guaranteed to be at least 
		// 1, so it's safe to add the last value.
		// This cannot be in the previous loop since
		// the shift shouldn't happen.
		tuple += TwoBit[*p];
		return 1;
	}

  int FromStringRL(Nucleotide *strPtr, TupleMetrics &tm) {

		//
		// Tuples are created with the right-most character
		// in the most significant bit in the sequence.
		//
		DNASequence tmpSeq;
		tmpSeq.seq = strPtr;
		tmpSeq.length = tm.tupleSize;
		if (!OnlyACTG(tmpSeq))
			return 0;

		if (tm.tupleSize == 0)
			return 1;

		tuple = 0;
		Nucleotide *p;
		for (p = strPtr + tm.tupleSize - 1; p > strPtr; p--) {
			tuple += TwoBit[*p];
			tuple <<=2;
		}
		//
		// The tuple size is guaranteed to be at least 
		// 1, so it's safe to add the last value.
		// This cannot be in the previous loop since
		// the shift shouldn't happen.
		tuple += TwoBit[*p];
		return 1;
	}
  
  int ShiftAddRL(Nucleotide nuc, TupleMetrics &tm) {
    if (ThreeBit[nuc] > 3) {
      return 0;
    }
    else {
      tuple >>= 2;
      tuple += (TwoBit[nuc] << ((tm.tupleSize-1)*2));
      return 1;
    }
  }
  
	int MakeRC(DNATuple &dest, TupleMetrics &tm) {
		int i;
		ULong tempTuple = tuple;
		dest.tuple = 0;
		ULong mask = 0x3;
		if (tm.tupleSize == 0)
			return 0;
		for (i = 0; i < tm.tupleSize - 1; i++ ) {
			dest.tuple += (~tempTuple & mask);
			tempTuple >>= 2;
			dest.tuple <<= 2;
		}
		dest.tuple += (~tempTuple & mask);
		return 1;
	}

	string ToString(TupleMetrics &tm) {
		int i;
		string s;
		ULong tempTuple = tuple;
		for (i = 0;i < tm.tupleSize;i++) {
			s.insert(s.begin(), TwoBitToAscii[tempTuple & 3]);
			tempTuple = tempTuple >> 2;
		}
		return s;
	}

};


class CompareByTuple {
 public:
  bool operator()(const DNATuple &lhs, const DNATuple &rhs) const {
    return lhs.tuple < rhs.tuple;
  }
};


class CountedDNATuple : public DNATuple {
 public:
  int count;
};

class PositionDNATuple : public DNATuple {
 public:
	PositionDNATuple& operator=(const PositionDNATuple &rhs) {
		pos = rhs.pos;
		tuple = rhs.tuple;
		return *this;
	}

	int operator<(const PositionDNATuple & pTuple) const {
		if (tuple < pTuple.tuple) 
			return 1;
		else if (tuple == pTuple.tuple) 
			return pos < pTuple.pos;
		else 
			return 0;
	}

	
	int operator==(const PositionDNATuple &pTuple) const {
		return tuple == pTuple.tuple and pos == pTuple.pos;
	}

	int operator<(const DNATuple &pTuple) const {
		return (tuple < pTuple.tuple);
	}

	int operator==(const DNATuple &pTuple) const {
		return tuple == pTuple.tuple;
	}
	
	int operator!=(const DNATuple &pTuple) const {
		return tuple != pTuple.tuple;
	}

  PositionDNATuple() : DNATuple() {
		pos = -1;
  }

	PositionDNATuple(PositionDNATuple &tupleP, DNALength posP) {
		tuple = tupleP.tuple;
		pos   = posP;
	}

};

class OrderPositionDNATuplesByPosition {
 public:
	int operator()(const PositionDNATuple &lhs, const PositionDNATuple &rhs) const {
		return lhs.pos < rhs.pos;
	}
};

class OrderPositionDNATuplesByTuple {
 public:
	int operator()(const PositionDNATuple &lhs, const PositionDNATuple &rhs) const {
		return lhs.tuple < rhs.tuple;
	}
};


//
// I don't have time to make this use good design, so I'm copying and barely modifying.
//

template<typename Sequence> 
int SearchSequenceForTuple(Sequence &seq, TupleMetrics &tm, DNATuple &queryTuple) {
	DNALength p;
	PositionDNATuple tempTuple, upperTuple;
	
	p = 0;
	DNALength tupleTransEnd = seq.length - tm.tupleSize + 1;
	DNALength cur = 0;
	DNALength curValidEnd = 0;

	//
	// Construct the mask-off bit pair for the shifted tuple.
	//
	PositionDNATuple maskLeftTuple;
	maskLeftTuple.tuple = 3;
	maskLeftTuple.tuple = maskLeftTuple.tuple << 2*tm.tupleSize;
	maskLeftTuple.tuple = ~maskLeftTuple.tuple;
	PositionDNATuple testTuple;
	while (curValidEnd < seq.length) {
		//
		// Search for the next available window that can be translated into a tuple.
		//
		cur = curValidEnd;
		while(curValidEnd < seq.length and IsACTG[seq.seq[curValidEnd]]) {
			curValidEnd++;
		}
		if (curValidEnd - cur >= tm.tupleSize) {
			//
			// Found a span that does not have N's in it, 
			//
			assert (tempTuple.FromStringRL(&(seq.seq[cur]), tm) == 1);
			p = cur;
			if (tempTuple.tuple == queryTuple.tuple) {
				return 1;
			}
			for (p++; p < curValidEnd - tm.tupleSize + 1; p++) {
				tempTuple.tuple >>=2;
				//				tempTuple.tuple &= maskLeftTuple.tuple;
				upperTuple.tuple = TwoBit[seq.seq[p+tm.tupleSize-1]];
				upperTuple.tuple = upperTuple.tuple << (2 * (tm.tupleSize-1));
				tempTuple.tuple += upperTuple.tuple;
				if (tempTuple.tuple == queryTuple.tuple) {
					return 1;
				}
			}
		}
		else {
			++curValidEnd;
		}
	}
}


template<typename Sequence>
int SequenceToTupleList(Sequence &seq, TupleMetrics &tm, TupleList<DNATuple> &tupleList) {
	DNALength p;
	PositionDNATuple tempTuple, upperTuple;
	
	p = 0;
	DNALength tupleTransEnd = seq.length - tm.tupleSize + 1;
	DNALength cur = 0;
	DNALength curValidEnd = 0;

	//
	// Construct the mask-off bit pair for the shifted tuple.
	//
	PositionDNATuple maskLeftTuple;
	maskLeftTuple.tuple = 3;
	maskLeftTuple.tuple = maskLeftTuple.tuple << 2*tm.tupleSize;
	maskLeftTuple.tuple = ~maskLeftTuple.tuple;
	PositionDNATuple testTuple;
	while (curValidEnd < seq.length) {
		//
		// Search for the next available window that can be translated into a tuple.
		//
		cur = curValidEnd;
		while(curValidEnd < seq.length and IsACTG[seq.seq[curValidEnd]]) {
			curValidEnd++;
		}
		if (curValidEnd - cur >= tm.tupleSize) {
			//
			// Found a span that does not have N's in it, 
			//
			assert (tempTuple.FromStringRL(&(seq.seq[cur]), tm) == 1);
			p = cur;
			tupleList.Append(tempTuple);
			for (p++; p < curValidEnd - tm.tupleSize + 1; p++) {
				tempTuple.tuple >>=2;
				//				tempTuple.tuple &= maskLeftTuple.tuple;
				upperTuple.tuple = TwoBit[seq.seq[p+tm.tupleSize-1]];
				upperTuple.tuple = upperTuple.tuple << (2 * (tm.tupleSize-1));
				tempTuple.tuple += upperTuple.tuple;
				//testTuple.FromStringRL(&seq.seq[p], tm);
				//assert(testTuple.tuple == tempTuple.tuple);
				tupleList.Append(tempTuple);
			}
		}
		else {
			++curValidEnd;
		}
	}
	return tupleList.size();
}


template<typename Sequence>
int SequenceToTupleList(Sequence &seq, TupleMetrics &tm, TupleList<PositionDNATuple> &tupleList) {
	DNALength p;
	PositionDNATuple tempTuple, upperTuple;
	
	p = 0;
	DNALength tupleTransEnd = seq.length - tm.tupleSize + 1;
	DNALength cur = 0;
	DNALength curValidEnd = 0;

	//
	// Construct the mask-off bit pair for the shifted tuple.
	//
	PositionDNATuple maskLeftTuple;
	maskLeftTuple.tuple = 3;
	maskLeftTuple.tuple = maskLeftTuple.tuple << 2*tm.tupleSize;
	maskLeftTuple.tuple = ~maskLeftTuple.tuple;
	PositionDNATuple testTuple;
	while (curValidEnd < seq.length) {
		//
		// Search for the next available window that can be translated into a tuple.
		//
		cur = curValidEnd;
		while(curValidEnd < seq.length and IsACTG[seq.seq[curValidEnd]]) {
			curValidEnd++;
		}
		if (curValidEnd - cur >= tm.tupleSize) {
			//
			// Found a span that does not have N's in it, 
			//
			assert (tempTuple.FromStringRL(&(seq.seq[cur]), tm) == 1);
			p = cur;
			tempTuple.pos = p;
			tupleList.Append(tempTuple);
			for (p++; p < curValidEnd - tm.tupleSize + 1; p++) {
				tempTuple.tuple >>=2;
				//				tempTuple.tuple &= maskLeftTuple.tuple;
				upperTuple.tuple = TwoBit[seq.seq[p+tm.tupleSize-1]];
				upperTuple.tuple = upperTuple.tuple << (2 * (tm.tupleSize-1));
				tempTuple.tuple += upperTuple.tuple;
				tempTuple.pos = p;
				//testTuple.FromStringRL(&seq.seq[p], tm);
				//assert(testTuple.tuple == tempTuple.tuple);
				tupleList.Append(tempTuple);
			}
		}
		else {
			++curValidEnd;
		}
	}
	return tupleList.size();
}



#endif
