#ifndef LIS_PVALUE_H_
#define LIS_PVALUE_H_

#include "../../datastructures/anchoring/MatchPos.h"
#include "../../datastructures/tuplelists/TupleCountTable.h"
#include "ScoreAnchors.h"

template<typename T_MatchPos>
void StoreNonOverlappingIndices(vector<T_MatchPos> &lis, vector<T_MatchPos> &noOvpLis) {

	unsigned int i;

	//
	// Greedily add lis matches according to weight.  A match may be added
	// as long as it does not overlap with any other matches.
	//

	// do nothing on empty lists
	if (lis.empty()) {
		return;
	}

	//
	// First build a list of matches sorted by weight.

	SortMatchPosListByWeight(lis);


	//
	// The first match is guaranteed to not overlap.
	noOvpLis.push_back(lis[0]);


	//
	// Nothing is overlapping, and everything is sorted when there is 
	// just one value.
	if (lis.size() == 1)
		return;

	//
	// Next, add matches as long as they do not overlap.
	for (i = 1; i < lis.size(); i++ ){

		VectorIndex j;
		int lts = lis[i].t;
		int lte = lis[i].t + lis[i].GetLength();
		int lqs = lis[i].q;
		int lqe = lis[i].q + lis[i].GetLength();

		int ovpFound = 0;
		for (j =0; j < noOvpLis.size(); j++ ){
			int tIntvStart = noOvpLis[j].t;
			int tIntvEnd   = noOvpLis[j].t + noOvpLis[j].GetLength();
			int qIntvStart = noOvpLis[j].q;
			int qIntvEnd   = noOvpLis[j].q + noOvpLis[j].GetLength();
			if ((lts >= tIntvStart and lts < tIntvEnd) or 
					(lte >  tIntvStart and lte <= tIntvEnd) or
					(lqs >= qIntvStart and lqs < qIntvEnd) or
					(lqe >  qIntvStart and lqe <= qIntvEnd)) {
				ovpFound = 1;
				break;
			}
		}
		if (!ovpFound) {
			noOvpLis.push_back(lis[i]);
		}
	}
	
	//
	// Now, the matches are found in order of size, but they need to
	// be stored in order of text.
	//
	SortMatchPosList(noOvpLis);

	//
	// The match pos list was sorted in order of weight. 
	// Just in case it causes problems down the line, re-sort it 
	// according to query pos.
	//
	lis = noOvpLis;
	SortMatchPosList(lis);
	
}


template<typename T_TextSequence, typename TSequence, typename T_MatchPos, typename T_Tuple>
	float ComputeLISPValue(vector<T_MatchPos> &lis, T_TextSequence &text, TSequence &read,
												 TupleMetrics &tm, TupleCountTable<T_TextSequence, T_Tuple> &ct,
                         int &lisNBases, int &lisSize ) {
	//
	// First, find a subset of the lis that has non-overlapping matches.
	//


  int i;
  lisNBases = 0;
  for (i = 0; i < lis.size(); i++) {
    lisNBases += lis[i].l;
  }
  lisSize = lis.size();

	float neginf = -1.0/0.0;
	float inf = 1.0/0.0;

	if (lis.size() == 1) {
		//
		// When just one LIS is found, the pvalue of the match is just the
		// pvalue of getting a single hit with one read, which we estimate
		// as the pvalue of at least one hit.
		// 
		float matchProb;
		//
		// Weight the single match based on how frequently it appears in the genome.

		if (POneOrMoreMatches<T_TextSequence>(text, lis[0].t, lis[0].GetLength(), tm, ct, matchProb)) {		
			return matchProb;
    }
    else {
			// return non-significant match value.
		  return 1;
    }
	}
	else {
		//
		// There is more than one non overlapping match.  This evokes a totally
		// different probability metric on LIS values.  Rather than computing the probability 
		// of a single match, compute the probability of seeing several matches in a row.
		//
		VectorIndex i;
		float pChain = 0; // In other words, pChain = log(1)
		// 
		// prob of chain starts out as prob of first match.
		//
		if (POneOrMoreMatches(text, lis[0].t, lis[0].GetLength(), tm, ct, pChain)) {
			//  assert(pChain != neginf);
			//
			// Next, for each match, compute the expected probability of having to wait 
			// at most tGap time for a match.  Since this has screened for anchors
			// that are not overlapping, consider the events to be independent, and therefore
			// the probabilities multiply.
			//
			for (i = 1; i < lis.size(); i++ ){ 

				int tupleMatchCount;
				//
				// Find out how many times this word appears in the genome.
				//
				int   approxNumMatches;
				//
				// Assume all hits are uniformly distributed across the
				// genome.  qLambda is the frequency of seeing hit i in the
				// genome.
				//
				float qLambda;
				int   tGap;
				int qGap;
				//				qLambda          = lis[i].GetMultiplicity() / (1.0*text.length);
				qLambda          = 1 / (1.0*text.length);
			
				// Now compute the probability of the chain.  Since the
				// matches are uniformly distributed across the genome, the
				// waiting time between matches is exponentially distributed.

				pChain = pChain + log(qLambda);
			}
			return pChain;
		}
		else {
			return 1;
		}
  }
}


#endif
