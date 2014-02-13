#ifndef TUPLE_METRICS_H_
#define TUPLE_METRICS_H_

#include "../Types.h"
#include "TupleMask.h"

class TupleMetrics {
 public:
	int tupleSize;
	ULong tupleMask;

  TupleMetrics() {
    tupleSize = tupleMask = 0;
  }

	void InitializeMask() {
    tupleMask = TupleMask[tupleSize];
	}

	void Initialize(int pTupleSize) {
		tupleSize = pTupleSize;
		InitializeMask();
	}

};


#endif
