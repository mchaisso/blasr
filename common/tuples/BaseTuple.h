#ifndef BASE_TUPLE_H_
#define BASE_TUPLE_H_

#include "../Types.h"
#include "TupleMetrics.h"

class BaseTuple {
 public:
  ULong  tuple;
	BaseTuple() {
		tuple = 0;
	}
  ULong  HashPowerOfFour(int nBases, TupleMetrics &tm) {
    //
    // When the hash can fit inside the entire tuple, just return the
    // tuple.
    //
    if (tm.tupleSize > nBases) {
      return tuple;
    }
    else {
      return ((tuple & TupleMask[nBases]) + (tuple % 1063)) % (1 << (nBases*2));
    }
  }
  /*  
  BaseTuple &operator=(const BaseTuple &rhs) {
		tuple = rhs.tuple;
		return *this;
  };
  */

  int operator<(const BaseTuple &rhs) const {
		return tuple < rhs.tuple;
  }

  int operator==(const BaseTuple &rhs) const {
		return tuple == rhs.tuple;
  }

	int operator!= (const BaseTuple &rhs) const {
		return tuple != rhs.tuple;
	}

	BaseTuple ShiftLeft(TupleMetrics &tm, ULong shift=1L) {
		tuple = tuple << shift;
		tuple = tuple & tm.tupleMask;
		return *this;
	}

	BaseTuple ShiftRight(ULong shift=1L) {
		tuple = tuple >> shift;
		return *this;
	}
	
	BaseTuple Append(ULong val, TupleMetrics &tm, ULong nBits) {
		tuple = tuple << nBits;
		tuple = tuple & tm.tupleMask;
		tuple = tuple + val;
		return *this;
	}

	long ToLongIndex() {
		long tempTuple = tuple;
		return tempTuple;
	}
};


#endif
