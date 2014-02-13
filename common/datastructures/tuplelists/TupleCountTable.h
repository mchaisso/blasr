#ifndef TUPLE_COUNT_TABLE_H_
#define TUPLE_COUNT_TABLE_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <assert.h>
#include "../../tuples/TupleMetrics.h"
using namespace std;

template<typename TSequence, typename TTuple>
class TupleCountTable {
 public:
	int *countTable;
	int countTableLength;
	int nTuples;
	TupleMetrics tm;
	bool deleteStructures;
	void InitCountTable(TupleMetrics &ptm) {
		tm = ptm;
		tm.InitializeMask();
		assert(tm.tupleSize > 0);
		// create the mask just in case the ptm is not initialized properly.
		countTableLength = 4;
		countTableLength = countTableLength << ((tm.tupleSize - 1)*2);
		InitCountTable(countTableLength);
	}

	TupleCountTable() {
		countTable = NULL;
		countTableLength = 0;
		nTuples = 0;
		deleteStructures = true;
	}

	~TupleCountTable() {
		if (deleteStructures == false) {
			//
			// Do not delete this if it is referencing another structure/
			//
			return;
		}
			
		if (countTable != NULL) {
			delete [] countTable;
		}
	}
	void InitCountTable(int p_countTableLength ){ 
		countTableLength = p_countTableLength;
		assert(countTableLength > 0);
		countTable = new int[countTableLength];
		fill(&countTable[0], &countTable[countTableLength], 0);
		nTuples = 0;
	}

	void IncrementCount(TTuple &tuple) {
		long tupleIndex = tuple.ToLongIndex();
		assert(tupleIndex < countTableLength);
		countTable[tupleIndex]++;
		++nTuples;
	}

	void AddSequenceTupleCountsLR(TSequence &seq) {
		VectorIndex i;
		TTuple tuple;
		if (seq.length>= tm.tupleSize) {
			for (i = 0; i < seq.length - tm.tupleSize + 1; i++ ){ 
				if (tuple.FromStringLR(&seq.seq[i], tm)) {
					IncrementCount(tuple);
				}
			}
		}
	}

	void Write(ofstream &out) {
		out.write((char*) &countTableLength, sizeof(int));
		out.write((char*) &nTuples, sizeof(int));
    out.write((char*) &tm.tupleSize, sizeof(int));
		out.write((char*) countTable, sizeof(int) * countTableLength);
	}
	
	void Read(ifstream &in) {
		in.read((char*) &countTableLength, sizeof(int));
		in.read((char*) &nTuples, sizeof(int));
    in.read((char*) &tm.tupleSize, sizeof(int));
		tm.InitializeMask();
    countTable = new int[countTableLength];
		in.read((char*) countTable, sizeof(int) * countTableLength);
	}
};


#endif
