#ifndef TUPLES_DNA_TUPLE_LIST_H_
#define TUPLES_DNA_TUPLE_LIST_H_

#include <vector>
#include "DNATuple.h"





template<typename Sequence>
	DNALength StoreTuplePosList(Sequence seq, TupleMetrics &tm, vector<PositionDNATuple> &tupleList) {
	//
	// Do this faster later on with a suffix tree -- faster than n log n construction time.
	// 
	DNALength s;
	
	PositionDNATuple tempTuple;
	for (s = 0; s < seq.length - tm.tupleSize + 1; s++) {
		if (tempTuple.FromStringLR(&(seq.seq[s]), tm)) {
			tempTuple.pos = s;
			tupleList.push_back(tempTuple);
		}
	}

	std::sort(tupleList.begin(), tupleList.end());
	
		// 
		// Be nice and leave the pos list in ascending sorted order,
		// even though the top of this function does not specify it.
		//
	
	return tupleList.size();
}


void WriteTuplePosList(vector<PositionDNATuple> &tupleList, int tupleSize, ofstream &out) {
	DNALength tupleListLength = tupleList.size();
	out.write((char*) &tupleSize, sizeof(tupleSize));
	out.write((char*) &tupleListLength, sizeof(DNALength));
	out.write((char*) &tupleList[0], sizeof(PositionDNATuple) * tupleList.size());
}

void ReadTuplePosList(ifstream &in, vector<PositionDNATuple> &tupleList, int &tupleSize) {
	DNALength tupleListLength;
	in.read((char*) &tupleSize, sizeof(int));
	in.read((char*) &tupleListLength, sizeof(DNALength));
	tupleList.resize(tupleListLength);
	in.read((char*) &tupleList[0], sizeof(PositionDNATuple) * tupleListLength);
}


#endif
