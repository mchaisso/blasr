#ifndef TUPLE_H_
#define TUPLE_H_

#include <vector>
#include <utility>

#include "Types.h"
#include "tuples/TupleMetrics.h"

template<typename Sequence, typename Tuple>
  BuildTupleList(Sequence seq, vector<Tuple> &tupleList) {
		BuildTupleList(seq, 0, seq.length, tupleList);
	}

template<typename Sequence, typename Tuple>
	BuildTupleList(Sequence seq, int seqStart, int length, vector<Tuple> &tupleList) {
		
	}

template<typename Sequence, typename CountedTuple>
	BuildCountedTupleList(Sequence seq, int seqStart, int length, TupleMetrics &tm, vector<CountedTuple> &tupleList) {
		int s;
		CountedTuple tuple;
		for (s = 0; s < seq.length - tm.tupleSize  + 1; s++) {
		}
	}


template<typename Sequence, typename Tuple>
	StoreUniqueTuplePosList(Sequence seq, TupleMetrics &tm, vector<int> &uniqueTuplePosList) {
		//
		// Do this faster later on with a suffix tree -- faster than n log n construction time.
		// 
		int s;
		vector<pair<Tuple, int> > tuples;
		Tuple tempTuple;
		for (s = 0; s < seq.length; s++) {
			tempTuple.FromStringLR(&(seq.seq[s]), tm);
			tuples.push_back(make_pair(tempTuple, s));
		}
		std::sort(tuples.begin(), tuples.end());
		int curUnique = 0, curPos = 0;

		//
		// Filter out the repetitive tuples.
		//

		while ( curPos < tuples.size() ) {
			int nextPos = curPos;

			while (nextPos < tuples.size() and tuples[nextPos] == tuples[curPos])
				nextPos++;
			if (nextPos - curPos == 1) {
				tuples[curUnique].first == tuples[curPos].first;
				uniqueTuplePosList.push_back(tuples[curUnique].second);
				++curUnique;
				++curPos;
			}
			else {
				curPos = nextPos;
			}
		}

		// 
		// Be nice and leave the pos list in ascending sorted order,
		// even though the top of this function does not specify it.
		//
		std::sort(uniqueTuplePosList.begin(), uniqueTuplePosList.end());
	}




	


				
			
		


#endif
