#ifndef MAP_BY_SUFFIX_ARRAY_H_
#define MAP_BY_SUFFIX_ARRAY_H_

#include "datastructures/suffixarray/SuffixArray.h"
#include "datastructures/anchoring/MatchPos.h"
#include "datastructures/anchoring/AnchorParameters.h"
#include "algorithms/alignment/SWAlign.h"
#include "algorithms/alignment/ScoreMatrices.h"

/*
 * Parameters:
 * Eventually this should be strongly typed, since this is specific to
 * suffix array searching on DNASequence read/genome types.
 * reference  - should be of type DNASequence
 * sa         - shuld be of type SuffixArray
 * read       - may be of any DNASequence subclass.
 * tupleSize  - The length of the keyword used to speed up searching.
 * Out:
 *   matchLow - The starting point in the suffix array for the LCP
 *              match for the read at pos p.
 *   matchHigh -The same array but for the upper bound.
 *   saMatchLength - The length of the lcp.
 */
template<typename T_SuffixArray, typename T_RefSequence, typename T_Sequence>
int LocateAnchorBoundsInSuffixArray(T_RefSequence &reference,
																		T_SuffixArray &sa,
																		T_Sequence &read, 
																		unsigned int minPrefixMatchLength,
																		vector<DNALength> &matchLow, 
																		vector<DNALength> &matchHigh,
																		vector<DNALength> &matchLength,
																		AnchorParameters &params) {

  //
  // Make sure there is enough of this read to map.  Since searches
  // are keyed off of 'minPrefixMatchLength' matches, don't search
  // anything shorter than that.
  //
  if (minPrefixMatchLength > 0 and 
      read.subreadEnd - read.subreadStart < minPrefixMatchLength) {
    return 0;
  }

	DNALength p, m;
	DNALength alignEnd;
  DNALength matchEnd = read.subreadEnd - minPrefixMatchLength + 1;
  DNALength numSearchedPositions = matchEnd - read.subreadStart;
  
	matchLength.resize(numSearchedPositions);
	matchLow.resize(numSearchedPositions);
	matchHigh.resize(numSearchedPositions);

  std::fill(matchLength.begin(), matchLength.end(), 0);
	std::fill(matchLow.begin(), matchLow.end(), 0);
	std::fill(matchHigh.begin(), matchHigh.end(), 0);
	vector<SAIndex> lowMatchBound, highMatchBound;	

	for (m = 0, p = read.subreadStart; p < matchEnd; p++, m++) {
		DNALength lcpLow, lcpHigh, lcpLength;
		lowMatchBound.clear(); highMatchBound.clear();
		lcpLow = 0;
		lcpHigh = 0;
		lcpLength = sa.StoreLCPBounds(reference.seq, reference.length, 
                                  &read.seq[p], matchEnd - p,
                                  params.useLookupTable,
                                  params.maxLCPLength,
                                  //
                                  // Store the positions in the SA
                                  // that are searched.
                                  //
                                  lowMatchBound, highMatchBound, 
                                  params.stopMappingOnceUnique);

    //
    // Possibly print the lcp bounds for debugging
    //
    if (params.lcpBoundsOutPtr != NULL) {
      int i;
      for (i = 0; i < lowMatchBound.size(); i++) {
        *params.lcpBoundsOutPtr << highMatchBound[i] - lowMatchBound[i];
        if (i < lowMatchBound.size() - 1) {
          *params.lcpBoundsOutPtr << " ";
        }  
      }
      *params.lcpBoundsOutPtr << endl;
    }

		//
		// Default to no match.
		//
		matchLow[m] = matchHigh[m] = matchLength[m] = 0;

    //
    // If anything was found in the suffix array:
    //
		if (lowMatchBound.size() > 0) {
      //
      // First expand the search bounds until at least one match is
      // found.
      //
      int lcpSearchLength = lowMatchBound.size();
      bool extendedForward = false;
      while (lcpSearchLength > 0 and 
             lowMatchBound[lcpSearchLength - 1] == 
             highMatchBound[lcpSearchLength - 1]) {
        lcpSearchLength--;
        lcpLength--;
      }
      matchLow[m]  = lowMatchBound[lcpSearchLength - 1];
      matchHigh[m] = highMatchBound[lcpSearchLength - 1];
      matchLength[m] = minPrefixMatchLength + lcpSearchLength - 1;

      //
      // Next, apply some heuristics to the anchor generation.
      //
      // 1.1 If the suffix array match is unique, try and extend that
      // match as long as possible to ease global chaining later on.  
      //
      // 1.2 If the suffix array match is unique, but cannot be
      // extended, it probably ends in an error.  Back the search up
      // by 1.
      //
      // 2.1 If the suffix array match is not unique, return the
      // default matches, or expand the search to include more
      // matches. 
      //

      //
      // Check to see if the match was unique.
      //

      if (matchLow[m] + 1 == matchHigh[m]) {
        //
        // If the match is unique, extend for as long as possible.
        //
        lcpLength = minPrefixMatchLength + lcpSearchLength - 1;
        long refPos    = sa.index[matchLow[m]] + lcpLength - 1;
        long queryPos  = p + lcpLength - 1;
        bool extensionWasPossible = false;

        while (refPos + 1 < reference.length and
               queryPos + 1 < read.length and
							 reference.seq[refPos + 1] != 'N' and 
               reference.seq[refPos + 1] == read.seq[queryPos + 1] and 
               (params.maxLCPLength == 0 or lcpLength < params.maxLCPLength)) {
          refPos++;
          queryPos++;
          lcpLength++;
          extensionWasPossible = true;
        }

        if (extensionWasPossible) {
          //
          // Was able to extend match far into the genome, store that.
          //
          matchLength[m] = lcpLength;
        }
        else if (extensionWasPossible == false) {
          //
          // No extension was possible, indicating that this match
          // ends at an error.  To be safe, expand search by up to
          // 1.
          //
          if (lcpSearchLength > 1) {
            lcpSearchLength = lcpSearchLength - 1;
          }
          matchLow[m]  = lowMatchBound[lcpSearchLength-1];
          matchHigh[m] = highMatchBound[lcpSearchLength-1];
          matchLength[m] = minPrefixMatchLength + lcpSearchLength - 1;
        }
      }
      else {
        //
        // The match is not unique.  Store a possibly expanded search.
        // 
        int numBacktrack = params.expand;
        if (lcpSearchLength > params.expand) {
          lcpSearchLength -= params.expand;
        }
        else {
          assert(lowMatchBound.size() > 0);
          lcpSearchLength = 1;
        }
          
        //
        // There are multiple matches for this position.
        //
        matchLow[m]    = lowMatchBound[lcpSearchLength - 1];
        matchHigh[m]   = highMatchBound[lcpSearchLength - 1];
        matchLength[m] = minPrefixMatchLength + lcpSearchLength - 1;
      }
    }
    else {
      //
      // The match is shorter than what the search is supposed to
      // expand to.  In order to avoid expanding to before the end
      // of the match list, do not set any match.
      //
      matchLow[m]    = 0;
      matchHigh[m]   = 0;
      matchLength[m] = 0;
		}
		if (params.advanceExactMatches) {
			//
			// Advance to past the mismatch.
			//
			p += max(1,(int)(lcpLength - params.advanceExactMatches));
			m += max(1,(int)(lcpLength - params.advanceExactMatches));
		}
	}
	return 1;
}

template<typename T_SuffixArray, typename T_RefSequence, typename T_Sequence, typename T_MatchPos>
int MapReadToGenome(T_RefSequence &reference,
										T_SuffixArray &sa,
										T_Sequence &read, 
										unsigned int minPrefixMatchLength,
										vector<T_MatchPos> &matchPosList,
										AnchorParameters &anchorParameters) {

	vector<DNALength> matchLow, matchHigh, matchLength;


  if (read.subreadEnd - read.subreadStart < anchorParameters.minMatchLength) {
    matchPosList.clear();
    return 0;
  }

	
	LocateAnchorBoundsInSuffixArray(reference, sa, read, minPrefixMatchLength, matchLow, matchHigh, matchLength, anchorParameters);

	//
	// Try evaluating some contexts.
	//
	DNALength pos;
	DNALength mappedLength = matchLow.size();
	assert(matchLow.size() == matchHigh.size());

	DNASequence evalQrySeq, evalRefSeq;
	vector<Arrow> pathMat;
	vector<int> scoreMat;
	Alignment alignment;
	
	
	//
	// Do some filtering on the matches looking for overlapping matches
	// if there are any.
	//
	if (anchorParameters.removeEncompassedMatches) {
		vector<bool> removed;
		removed.resize(read.length);
		std::fill(removed.begin(), removed.end(), false);
		int i;
		int nRemoved = 0;
		for (i = 0; i < read.length-1; i++) {
			if (matchLength[i] == matchLength[i+1]+1) {
				removed[i+1] = true;
			}
		}
		for (i = 1; i < matchLength.size(); i++) {
			if (removed[i]) {
				matchLength[i] = matchLow[i] = matchHigh[i] = 0;
			}
		}
	}
  //
  // Now add 
  // 
  DNALength endOfMapping;
  DNALength trim = max(anchorParameters.minMatchLength + 1, sa.lookupPrefixLength + 1);
  if (read.subreadEnd < trim) {
    endOfMapping = 0;
  }
  else {
    endOfMapping = read.subreadEnd - trim;
  }
	int totalMatched = 0;
 	for (pos = read.subreadStart; pos < endOfMapping; pos++) {	
    int matchIndex = pos - read.subreadStart;
    assert(matchIndex < matchHigh.size());
		if (matchHigh[matchIndex] - matchLow[matchIndex] <= anchorParameters.maxAnchorsPerPosition) {
			DNALength mp;
			for (mp = matchLow[matchIndex]; mp < matchHigh[matchIndex]; mp++ ) {
				if (matchLength[matchIndex] < anchorParameters.minMatchLength) {
					continue;
				}

				//
				// By default, add all anchors.
				//
				if (matchLength[matchIndex] + pos > read.length) {
					//
					// When doing branching, it's possible that a deletion
					// branch finds an anchor that goes past the end of a
					// read.  When that is the case, trim back the anchor
					// match since this confuses downstream assertions.
					//
					matchLength[matchIndex] = read.length - pos;
				}
        assert(sa.index[mp] + matchLength[matchIndex] <= reference.length);
				assert(reference.seq[sa.index[mp] + matchLength[matchIndex] - 1] != 'N');
				matchPosList.push_back(ChainedMatchPos(sa.index[mp], pos, matchLength[matchIndex]));
			}
		}
	}

	return matchPosList.size();

}		


										 
#endif
