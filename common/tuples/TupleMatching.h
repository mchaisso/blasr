#ifndef TUPLE_MATCHING_H_
#define TUPLE_MATCHING_H_

#include <vector>
#include <utility>
#include <iostream>

#include "TupleList.h"
#include "TupleMetrics.h"
#include "DNATuple.h"


template<typename Sequence, typename T_TupleList> 
	int SequenceToTupleList(Sequence &seq, TupleMetrics &tm, T_TupleList &tupleList) {
	int s;
	typename T_TupleList::Tuple tempTuple;
  if (seq.size() < tm.tupleSize) {
    return 1; 
  }
  
  // Otherwise, there is at least one tuple

  tupleList.Append(tempTuple);
  int res = 0;
	for (s = 0; s < seq.length - tm.tupleSize + 1; s++ ) {
    if ((res and (res = tempTuple.ShiftAddRL(seq.seq[s+tm.tupleSize-1], tm))) or
        (!res and (res = tempTuple.FromStringRL(&seq.seq[s], tm)))) {
      tempTuple.ShiftAddRL(seq.seq[s + tm.tupleSize - 1], tm);
      tempTuple.pos = s;
      tupleList.Append(tempTuple);
    }
	}
	return 1;
}
	

template<typename TSequence, typename TMatch, typename T_TupleList>
	int StoreMatchingPositions(TSequence &querySeq, TupleMetrics &tm, T_TupleList &targetTupleList, vector<TMatch> &matchSet, int maxMatches=0) {
	DNALength s;
//	TQueryTuple queryTuple;
	typename T_TupleList::Tuple queryTuple;
  queryTuple.pos = 0;
	if (querySeq.length >= tm.tupleSize) {
    int res = 0;
		for (s = 0; s < querySeq.length - tm.tupleSize + 1; s++) {
      if ((res and (res = queryTuple.ShiftAddRL(querySeq.seq[s+tm.tupleSize-1], tm))) or
          (!res and (res = queryTuple.FromStringRL(&querySeq.seq[s], tm)))) {
				int targetListIndex = 0;
        typename vector<typename T_TupleList::Tuple>::const_iterator curIt, endIt;
				targetTupleList.FindAll(queryTuple, curIt, endIt);
				if (maxMatches == 0 or endIt - curIt <= maxMatches) {
					for(; curIt != endIt; curIt++) {
						matchSet.push_back(TMatch(s, (*curIt).pos));
						++targetListIndex;
					}
				}
			}
		}
	}
	return matchSet.size();
}


template<typename Sequence, typename Tuple>
	int StoreUniqueTuplePosList(Sequence seq, TupleMetrics &tm, vector<int> &uniqueTuplePosList) {
		//
		// Do this faster later on with a suffix tree -- faster than n log n construction time.
		// 
		int s;
		vector<pair<Tuple, int> > tuples;
		Tuple tempTuple;
		for (s = 0; s < seq.length - tm.tupleSize + 1; s++) {
			tempTuple.FromStringRL(&(seq.seq[s]), tm);
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
		return uniqueTuplePosList.size();
	}

#endif
