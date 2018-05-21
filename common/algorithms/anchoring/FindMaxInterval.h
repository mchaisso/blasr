#ifndef MAX_INCREASING_INTERVAL
#define MAX_INCREASING_INTERVAL

#include <semaphore.h>
#include <fstream>
#include <iostream>
#include "GlobalChain.h"
#include "SDPGlobalChain.h"
#include "BasicEndpoint.h"
#include "IntervalSearchParameters.h"
#include "datastructures/anchoring/WeightedInterval.h"
#include "datastructures/anchoring/MatchPos.h"

template<typename T_Sequence, typename T_AnchorList>
class DefaultWeightFunction {
 public:
	float operator()(T_Sequence &text, T_Sequence &read, T_AnchorList matchPosList) {
		int i;
		float weight = 0;
		for (i = 0; i < matchPosList.size(); i++) {
			weight += matchPosList[i].weight();
		}
		return weight;
	}
};


template<typename T_Pos>
class MatchPosQueryOrderFunctor {
 public:
	int operator()(T_Pos &pos) {
		return pos.q;
	}
};

unsigned int NumRemainingBases(DNALength curPos, DNALength intervalLength) {
  if (curPos > intervalLength) {
    return 0;
  }
  else {
    return intervalLength - curPos;
  }
}


template<typename T_MatchList>
void RetainLongestContiguousMatch(T_MatchList &matchList, int maxGap) {
	int i, j;
	int maxSpanStart = 0;
	int maxSpanEnd   = 0;
	int maxSpanWeight = 0;
  //
  // First find the locations of the large gaps
  //
  int prevGapEnd = 0;
	int maxMatchStart = 0;
	int maxMatchEnd = 0;
	int maxMatchAnchorWeight = 0;
	int curMatchAnchorWeight = matchList[0].l;

  for (i = 0; i < matchList.size() - 1; i++) {
		if (matchList[i+1].t - matchList[i].t > maxGap or matchList[i+1].q - matchList[i].q > maxGap  ) {
			//
			// Gap between i and i+1
			//
			if (curMatchAnchorWeight > maxMatchAnchorWeight) {
				maxMatchStart = prevGapEnd;
				maxMatchEnd   = i;
				maxMatchAnchorWeight = curMatchAnchorWeight;
			}
			prevGapEnd = i+1;
			curMatchAnchorWeight = 0;
		}
		curMatchAnchorWeight += matchList[i+1].l;
	}
	if (curMatchAnchorWeight > maxMatchAnchorWeight) {
		maxMatchEnd = matchList.size();
		maxMatchStart = prevGapEnd;
	}
	
	for (i = 0; i < maxMatchEnd - maxMatchStart; i++) {
		matchList[i] = matchList[maxMatchStart+i];
	}

	matchList.resize(maxMatchEnd - maxMatchStart);
}

template<typename T_MatchList>
void ScreenSpuriousMatches(T_MatchList &matchList, 
													 int maxGap) {
	vector<bool> spurious(matchList.size(), false);
	int spurGap = maxGap / 4;
	int i, j;
  for (i = 1; i < matchList.size() - 1; i++) {
		//
		// If the distance between two anchors is sufficiently small, this may initiate a match
		//
		if ((matchList[i].t - matchList[i-1].t > spurGap and matchList[i+1].q - matchList[i].q > spurGap ) or 
				(matchList[i].q - matchList[i-1].q > spurGap and matchList[i+1].t - matchList[i].t > spurGap )) {
			spurious[i] = true;
		}
	}
	j = 0;
	for (i = 0; i < matchList.size(); i++) {
		if (spurious[i] == false) {
			matchList[j] = matchList[i];
			j++;
		}
	}

	matchList.resize(j);
}
template<typename T_MatchList>
void DivideIntoContiguousMatches(T_MatchList &matchList, 
																 int maxGap, 
																 int minLength,
																 vector<T_MatchList> &dividedLists) {
	ScreenSpuriousMatches(matchList, maxGap);
	int i, j;
	vector<int> matchStart, matchEnd;
	bool inMatch = false;
	int start=0, end;
  for (i = 0; i < matchList.size() - 1; i++) {
		//
		// If the distance between two anchors is sufficiently small, this may initiate a match
		//
		if (matchList[i+1].t - matchList[i].t < maxGap and
				matchList[i+1].q - matchList[i].q < maxGap  ) {
			if (inMatch == false) {
				start = i;
				inMatch = true;
			}
		}

		if (matchList[i+1].t - matchList[i].t > maxGap or
				matchList[i+1].q - matchList[i].q > maxGap) {
			if (inMatch == true) {
				end = i;
				if (end - start > minLength) {
					matchStart.push_back(start);
					matchEnd.push_back(end);
				}
			}
			inMatch=false;
		}
	}
	if (inMatch == true) {
		matchStart.push_back(start);
		matchEnd.push_back(matchList.size());
	}

	int c;
	for (i = 0; i < matchStart.size(); i++) {
		int l = matchEnd[i]-  matchStart[i];
		dividedLists.push_back(T_MatchList(l));
		for (j = 0; j < l; j++) {
			dividedLists[i][j] = matchList[matchStart[i] + j];
		}
	}
}



template<typename T_MatchList>
int MatchListWeight(T_MatchList &matchList) {
	int i;
	int totalWeight = 0;
	for (i = 0; i < matchList.size(); i++) {
		totalWeight+= matchList[i].l;
	}
	return totalWeight;
}



template<typename T_MatchList>
void PrintLIS(T_MatchList &matchList, ofstream &out) {
  int i;
	for (i = 0; i < matchList.size(); i++) {
		out << matchList[i].q << "\t" << matchList[i].t << "\t" << matchList[i].l << endl;
	}
  out << endl << endl;
}


template<typename T_MatchList,
         typename T_SequenceDB>
  void FilterMatchesAsLIMSTemplateSquare(T_MatchList &matches, 
                                         DNALength queryLength,
                                         DNALength limsTemplateLength,
                                         T_SequenceDB &seqDB) {
  int seqIndex;
  //
  // Make sure there is sequence coordinate information.
  //
  if (seqDB.nSeqPos == 0) {
    return;
  }
  int matchIndex = 0;
  for (seqIndex = 1; seqIndex < seqDB.nSeqPos; seqIndex++) {
    DNALength refLength = seqDB.seqStartPos[seqIndex] - seqDB.seqStartPos[seqIndex - 1];
    refLength = queryLength * 1.15 + limsTemplateLength; // account for indel error.
    //
    // Flag matches that are beyond the (rough) square with the length
    // of the query for removal.
    //
    while (matchIndex < matches.size() and 
           matches[matchIndex].t < seqDB.seqStartPos[seqIndex]) {
      if (matches[matchIndex].t > seqDB.seqStartPos[seqIndex-1] + refLength) {
        matches[matchIndex].l = 0;
      }
      matchIndex++;
    }
    int curMatchIndex = 0;
    matchIndex = 0;
    for (matchIndex = 0; matchIndex < matches.size(); matchIndex++) {
      if (matches[matchIndex].l != 0) {
        matches[curMatchIndex] = matches[matchIndex];
        curMatchIndex++;
      }
    }
    matches.resize(curMatchIndex);
  }
 }
         
 template<typename T_MatchList,
      	 typename T_SequenceBoundaryDB>
  void AdvanceIndexToPastInterval(T_MatchList &pos, DNALength nPos,
                                  DNALength intervalLength,
                                  T_SequenceBoundaryDB &SeqBoundary,
                                  DNALength startIndex, DNALength startIntervalBoundary,
                                  DNALength &index, DNALength &indexIntervalBoundary
                                  ) {

	 //
	 // Move the index forward until pos[index] is past 
	 // pos[startIndex]+length of read, or in a new contig.
	 //

  if (index >= pos.size()) {
    return;
  }
  indexIntervalBoundary = SeqBoundary(pos[index].t);
  DNALength boundaryIndex = SeqBoundary.GetIndex(pos[index].t);
  DNALength nextBoundary  = SeqBoundary.GetStartPos(boundaryIndex + 1);

  while (// index is not past the end of the genome
         index < nPos and 
         //
         // Stop when the index goes too far ahead.
         //
         pos[index].t - pos[startIndex].t < intervalLength and
         //
         // Still searching in the current contig.
         //
         indexIntervalBoundary == startIntervalBoundary) {
		index++;
		if (index < nPos) {
		  indexIntervalBoundary = SeqBoundary(pos[index].t);
		}
  }
}

template<typename T_MatchList>  
int RemoveZeroLengthAnchors(T_MatchList &matchList) {       
  int origSize = matchList.size();
  int cur = 0, m;
  for (m = 0; m < matchList.size(); m++) {
    if (matchList[m].l > 0) {
      matchList[cur] = matchList[m];
      cur++;
    }
  }
  matchList.resize(cur);
  return origSize - cur;
}

template<typename T_MatchList>
int RemoveOverlappingAnchors(T_MatchList &matchList) {       
  int cur = 0;
  int m;
  int n;
	if (matchList.size() == 0) {
		return 0;
	}

  for (m = matchList.size() - 1; m > 0; m--) {
    n = m - 1;
    //
    // Skip past repeats in the query.
    while (n > 0 and matchList[n].t == matchList[m].t) {
      n--;
    }

    bool mergeFound = false;
    int ni = n;
    while (mergeFound == false and 
           n > 0 and 
           matchList[n].t == matchList[ni].t) {
      if (matchList[n].q < matchList[m].q and
          matchList[n].t < matchList[m].t and
          matchList[n].l + matchList[n].q >= matchList[m].l + matchList[m].q and
          matchList[n].l + matchList[n].t >= matchList[m].l + matchList[m].t) {
        matchList[m].l = 0;
        mergeFound = true;
      }
      n--;
    }
  }
  int numRemoved = RemoveZeroLengthAnchors(matchList);
  return numRemoved;
}
template<typename T_MatchList>
int SumAnchors(T_MatchList &pos, int start, int end) {
	int sum = 0;
	int i;
	for (i = start; i < end; i++) {
		sum += pos[i].l;
	}
	return sum;
}

template<typename T_MatchList,
      	 typename T_SequenceBoundaryDB>
	void AdvanceOverlap(T_MatchList &pos, 
											VectorIndex nPos,
											// End search for intervals at boundary positions
											// stored in seqBoundaries
											DNALength intervalLength,  
											// How many sets to keep track of
											int minSize,
											
											T_SequenceBoundaryDB & ContigStartPos,
											IntervalSearchParameters &params,
											VectorIndex &cur, VectorIndex &next,
											DNALength &curBoundary, 
											DNALength &nextBoundary,
											vector<DNALength> &start,
											vector<DNALength> &end) {

	while (next < nPos) {
		//
		// Only need to find overlaps at the beginning and end of every contig.
		//
		curBoundary = ContigStartPos(pos[cur].t);

		
		AdvanceIndexToPastInterval(pos, nPos, intervalLength, ContigStartPos,
															 cur, curBoundary, next, nextBoundary);

		if (nextBoundary != curBoundary and nextBoundary - curBoundary < float(intervalLength)/2.0) {
			while (next < nPos and nextBoundary == ContigStartPos(pos[next+1].t)) {
				next++;
			}
			cur = next;
			next += 1;
			continue;
		}

		if (curBoundary == ContigStartPos(pos[next-1].t)) {
			if (next - cur > minSize) {
				start.push_back(cur);
				end.push_back(next);
			}
		}

		if (curBoundary == nextBoundary) {
			while (curBoundary  == nextBoundary) {
				next++;
				nextBoundary  = ContigStartPos(pos[next].t);
			}
			while (cur < nPos and pos[cur].t + intervalLength < pos[next-1].t) {
				cur++;
			}
			curBoundary = ContigStartPos(pos[cur].t);
			if (cur < nPos and curBoundary == ContigStartPos(pos[next-1].t )) {
				if (next - cur > minSize) {
					start.push_back(cur);
					end.push_back(next);
				}
			}
		}
		cur = next;
		next ++;
		nextBoundary = ContigStartPos(pos[next].t);
	}
}





template<typename T_MatchList,
      	 typename T_SequenceBoundaryDB>
	void StoreLargestIntervals(T_MatchList &pos, 
							   // End search for intervals at boundary positions
							   // stored in seqBoundaries
							   T_SequenceBoundaryDB & ContigStartPos,
										
							   //
							   // parameters
							   // How many values to search through for a max set.

							   DNALength intervalLength,  
							   // How many sets to keep track of
							   int minSize,
							   vector<DNALength> &start,
							   vector<DNALength> &end,
							   IntervalSearchParameters &params,
							   const char* name=NULL) {
	if (pos.size() == 0) {
		return;
	}
	if (pos.size() == 1) {
		start.push_back(0);
		end.push_back(0);
		return;
	}
		
  //
  // Search for clusters of intervals within the pos array within
  // pos[cur...next).  The value of 'next' should be the first anchor
  // outside the possible range to cluster, or the end of the anchor list.
	VectorIndex cur = 0;
	VectorIndex nPos = pos.size();
	VectorIndex endIndex = cur + 1;
	DNALength curBoundary = 0, endIndexBoundary = 0;
	
	DNALength selfStart = 0;
	DNALength selfEnd   = 0;
	  
	if (params.noSelf == true) {
	  ContigStartPos.GetBoundaries(name, selfStart, selfEnd);
	}

	if (params.overlap == true) {
		AdvanceOverlap(pos, nPos, 
					   intervalLength, minSize, 
					   ContigStartPos, 
					   params, cur, endIndex, curBoundary, endIndexBoundary, start, end);
	}
	else { 
		curBoundary = ContigStartPos(pos[cur].t);
		if (endIndex < pos.size()) {
			endIndexBoundary = ContigStartPos(pos[endIndex].t);
		}
		else {
			endIndexBoundary = ContigStartPos(pos[cur].t);
		}

		//
		// Advance endIndex until the anchor is outside the interval that
		// statrts at 'cur', and is inside the same contig that the anchor
		// at cur is in.
		//

		DNALength curIntervalLength = NumRemainingBases(pos[cur].q, intervalLength);
		if (params.noSelf) { 
		  DNALength beforeSelfAdvance=cur;
		  while ( cur < nPos and pos[cur].t >= selfStart and pos[cur].t <= selfEnd) {
			cur+=1;
		  }
		  endIndex = cur+1;
		  curBoundary = ContigStartPos(pos[cur].t);
		  endIndexBoundary = ContigStartPos(pos[endIndex].t);
		  
		}
		AdvanceIndexToPastInterval(pos, nPos, intervalLength, ContigStartPos,
								   cur, curBoundary, endIndex, endIndexBoundary);



		DNALength prevStart = cur, prevEndIndex = endIndex ;
		int prevSize = endIndex - cur;
		DNALength maxStart = cur, maxEndIndex = endIndex;
		int maxSize = SumAnchors(pos, cur, endIndex);
		int curSize = maxSize;


		  
		
		if (curSize > minSize) {
			start.push_back(cur);
			end.push_back(endIndex);
		}
		while ( cur < nPos ) {
			// 
			// This interval overlaps with a possible max start
			//

			if (pos[cur].t >= pos[maxStart].t and 
					maxEndIndex > 0 and 
					pos[cur].t < pos[maxEndIndex].t and 
					curBoundary == endIndexBoundary) {

				if (curSize > maxSize) {
					maxSize = curSize;
					maxStart = cur;
					maxEndIndex   = endIndex;
				}
			}
			else {
				if (maxSize > minSize) {
					start.push_back(maxStart);
					end.push_back(maxEndIndex);
				}
				maxStart = cur;
				maxEndIndex   = endIndex;
				maxSize  = curSize;
			}
		
			//
			// Done scoring current interval.  At this point the range
			// pos[cur...endIndex) has been searched for a cluster of matches.
			// Find a new range that will possibly yield a new
			// maximum interval.  
  
			// 
			// If the endIndex anchor is not within the same contig as the current
			// anchor, advance cur to endIndex since it is impossible to find
			// any more intervals in the current pos.  Because curSize is a 
			// sliding window, the size needs to get reset.
			//
			if (curBoundary != endIndexBoundary) {
				cur = endIndex;
				curBoundary = endIndexBoundary;

				//
				// Start the search for the first interval in the endIndex contig
				// just after the current position.
				//
				if (endIndex < nPos) {
					endIndex++;
					if (params.noSelf) { 
					  DNALength beforeSelfAdvance=cur;
					  while ( cur < nPos and pos[cur].t >= selfStart and pos[cur].t <= selfEnd) {
						cur+=1;
					  }
					  endIndex = cur+1;
					  curBoundary = ContigStartPos(pos[cur].t);
					  endIndexBoundary = ContigStartPos(pos[endIndex].t);
					}

					AdvanceIndexToPastInterval(pos, nPos, intervalLength, ContigStartPos,
																		 cur, curBoundary, endIndex, endIndexBoundary);
					maxSize = 0;
					maxStart = cur;
					maxEndIndex = endIndex;
				}
				// Reset the sliding window size.
				maxSize = curSize = SumAnchors(pos, cur, endIndex);

			}
			else {
				
				//
				// The endIndex anchor is in the same contig as the current anchor
				//

				//
				// The interval [cur,endIndex) has been scored.  Move the cur
				// index forward until the interval contains the endIndex match,
				// and score that.

				// In ascii art:
				// genome  |---+----+---+-+---+--------+-------------+-----+---|
				//            cur  cur+1                             endIndex endIndex+1
				// read intv   ==================== 
				//
				// Movest to:
				// genome  |---+----+---+-+---+--------+-------------+-----+---|
				//            cur  cur+1                             endIndex endIndex+1
				// read intv                           ==================== 

				//
				// Because curSize is a sliding window, remove the anchors
				// that are skipped from the window size.
				//

				while (cur < endIndex and 
							 pos[cur].t + intervalLength < pos[endIndex].t) {
					curSize -= pos[cur].l;
					cur++;
				}
					
				//
				// At last anchor. Outside this loop will possibly add the
				// current interval to the set of candidates if it is large
				// enough. 
				//
				if (cur >= nPos)
					break;
					
				//
				// Advance the endIndex to outside this interval.
				//
				curSize += pos[endIndex].l;
				endIndex++;
			}

			//			cout << "scan: " << endIndex << " " << curSize << endl;
			if (endIndex > nPos) {
				//
				// Searched last interval, done.
				//
				break;
			}
	
			//
			//  When searching multiple contigs, it is important to know the
			//  boundary of the contig that this anchor is in so that clusters
			//  do not span multiple contigs.  Find the (right hand side)
			//  boundary of the current contig.
			//

			curBoundary  = ContigStartPos(pos[cur].t);
			endIndexBoundary = ContigStartPos(pos[endIndex].t);  
    
			//
			// Previously tried to advance half.  This is being removed since
			// proper heuristics are making it not necessary to use.
			//
		}

		if (curSize > minSize) {
			start.push_back(maxStart);
			end.push_back(maxEndIndex);
		}
	}
}

template <typename T_MatchList>
void SwapXY(T_MatchList &pos, int start, int end) {
	int i;
	int q;
	for (i = start; i < end;i++){
		q=pos[i].q;
		pos[i].q=pos[i].t;
		pos[i].t=q;
	}
}

template<typename T_MatchList>
	void StoreLIS(T_MatchList &pos, vector<VectorIndex> lisIndices, DNALength cur, T_MatchList &lis) {
	int i;
	// Maybe this should become a function?
	for (i = 0; i < lisIndices.size(); i++) {	
		lis.push_back(pos[lisIndices[i]+cur]); 
	}
}


template<typename T_MatchList,
	       typename T_PValueFunction, 
	       typename T_WeightFunction,
      	 typename T_SequenceBoundaryDB,
         typename T_ReferenceSequence,
	       typename T_WeightedIntervalSet,
	       typename T_Chained_Anchor,
	       typename T_Sequence>
	int FindMaxIncreasingInterval(//
																// Input
																// readDir is used to indicate if the interval that is being stored is in the forward
																// or reverse strand.  This is important later when refining alignments so that the
																// correct sequene is aligned back to the reference.
																// 
																int readDir, 
																T_MatchList &pos, 

																//
																// parameters
																// How many values to search through for a max set.

																DNALength intervalLength,  
																// How many sets to keep track of
																VectorIndex nBest, 

																// End search for intervals at boundary positions
																// stored in seqBoundaries
																T_SequenceBoundaryDB & ContigStartPos,
																
																// First rand intervals by their p-value
																T_PValueFunction &MatchPValueFunction,  

																// When ranking intervals, sum over weights determined by MatchWeightFunction
																T_WeightFunction &MatchWeightFunction,  

																//
																// Output.
																// The increasing interval coordinates, 
																// in order by queue weight.
																T_WeightedIntervalSet &intervalQueue,
																T_ReferenceSequence &reference, T_Sequence &query,
																IntervalSearchParameters &params,
																vector<BasicEndpoint<T_Chained_Anchor> > *chainEndpointBuffer,
																vector<Fragment> *fragmentBuffer,
                                const char *titlePtr=NULL
																) {

	int nRectangles;
	multiset<WeightedInterval<T_Chained_Anchor>, CompareWeightedIntervalByPValue<T_Chained_Anchor> > sdpiq;
	VectorIndex cur = 0;
	VectorIndex nPos = pos.size();
	vector<VectorIndex> lisIndices;
	//
	// Initialize the first interval.
	//
	if (pos.size() == 0) {
		return 0;
	}

	int lisSize;
	float lisWeight;
	float lisPValue;
	T_MatchList lis;
	float neginf = -1.0/0.0;	
  int noOvpLisSize = 0;
  int noOvpLisNBases = 0;
  
  // cut from here 1

  //
  // Search for clusters of intervals within the pos array within
  // pos[cur...end).  The value of 'end' should be the first anchor
  // outside the possible range to cluster, or the end of the anchor list.

	VectorIndex endIndex = cur + 1;
	DNALength curBoundary = 0, endIndexBoundary = 0;
	vector<UInt> scores, prevOpt;

	vector<DNALength> start, end;
	
	StoreLargestIntervals(pos, ContigStartPos, intervalLength, params.minInterval, start, end, params, titlePtr);
	VectorIndex i;
	VectorIndex posi;
	int maxLISSize = 0;
	//	cerr << "Stored " << start.size() << " largest intervals." << endl;
	for (posi = 0; posi < start.size(); posi++) {
		
		lis.clear();
		lisIndices.clear();
		cur = start[posi];
		endIndex = end[posi];
		/*		cerr << "Range t: " << pos[start[posi]].t << " - " << pos[end[posi]].t << "  q: "
					<< pos[start[posi]].q << " - " << pos[end[posi]].q << " " << end[posi] - start[posi] << endl;*/
		if (endIndex - cur == 1) {
      //
      // Just one match in this interval, don't invoke call to global chain since it is given.
      //
			lisSize = 0;
			lisIndices.push_back(0);
		}
		else {
      //
      // Find the largest set of increasing intervals that do not overlap.
      //
			string name;

			if (params.globalChainType == 0) {
				lisSize = GlobalChain<T_Chained_Anchor, BasicEndpoint<T_Chained_Anchor> >(pos, cur, endIndex, 
																																								lisIndices, chainEndpointBuffer);
			}
			else {
        //
        //  A different call that allows for indel penalties.
        //
				lisSize = SDPGlobalChain(&pos[cur], endIndex-cur, lisIndices, params.minMatch, *fragmentBuffer);
			}
			
		}

		//
		// Store the subset of pos in the longest increasing subsequence in lis
		//
		vector<T_MatchList> contiguousMatches;
		StoreLIS(pos, lisIndices, cur, lis);
		if (params.maxAnchorGap > 0) {
			//
			// This can divide a single hit into multiple smaller ones.
			//

			DivideIntoContiguousMatches(lis,
																	params.maxAnchorGap,
																	100,
																	contiguousMatches);

		}
		else {
			contiguousMatches.resize(1);
			contiguousMatches[0] = lis;
			//			cerr << "checking lis of size: " << lis.size() <<endl;
		}

		int cm;
		for (cm = 0; cm < contiguousMatches.size(); cm++) {
			
			
			// 
			// Compute pvalue of this match.
			//
			T_MatchList *lisPtr = &contiguousMatches[cm];
			lisSize = MatchListWeight(*lisPtr);

			if (lisPtr->size() > 0) {
				lisPValue = MatchPValueFunction.ComputePValue(*lisPtr,
																											noOvpLisNBases, 
																											noOvpLisSize);
			}
			else {
				lisPValue = 0;
			}

			if (lisSize > maxLISSize) {
				maxLISSize  = lisSize;
			}

			//
			// Insert the interval into the interval queue maintaining only the 
			// top 'nBest' intervals. 
			//
			noOvpLisSize = (*lisPtr).size();
			int lisIndex = 0;
			for (lisIndex = 0; lisIndex < (*lisPtr).size(); lisIndex++) {
				noOvpLisNBases+=(*lisPtr)[lisIndex].l;
			}

			typename WeightedIntervalSet<T_Chained_Anchor>::iterator lastIt = intervalQueue.begin();
			MatchWeight lisWeight = MatchWeightFunction(*lisPtr);
			VectorIndex lisEnd    = lisPtr->size() - 1;
			//			cerr << "nbases: " << noOvpLisNBases << " pvalue: " << lisPValue << " max: " << params.maxPValue << endl;
			
			if (lisPValue < params.maxPValue and lisSize > 0 and noOvpLisNBases > params.minInterval  ) {
				
				WeightedInterval<T_Chained_Anchor> weightedInterval(lisWeight, noOvpLisSize, noOvpLisNBases, 
																					(*lisPtr)[0].t, (*lisPtr)[lisEnd].t + (*lisPtr)[lisEnd].GetLength(), 
																					readDir, lisPValue, 
																					(*lisPtr)[0].q, (*lisPtr)[lisEnd].q + (*lisPtr)[lisEnd].GetLength(), query.length,
																					*lisPtr);
				intervalQueue.insert(weightedInterval);
			}
		}
	}
	return maxLISSize;
}


#endif
