#ifndef MAX_INCREASING_INTERVAL
#define MAX_INCREASING_INTERVAL

#include <semaphore.h>
#include <fstream>
#include <iostream>
#include "GlobalChain.h"
#include "SDPGlobalChain.h"
#include "BasicEndpoint.h"
#include "datastructures/anchoring/WeightedInterval.h"
#include "datastructures/anchoring/MatchPos.h"
#include "datastructures/anchoring/ClusterList.h"
#include "statistics/VarianceAccumulator.h"

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

class IntervalSearchParameters {
 public:
	int   globalChainType;
	float maxPValue;
	bool  overlap;
	int   minMatch;
	int   minInterval;
	IntervalSearchParameters() {
		globalChainType = 0;
		maxPValue       = log(0.1);
		overlap         = false;
		minMatch        = 0;
		minInterval     = 0;
	}
};


template<typename T_MatchList>
void PrintLIS(T_MatchList &matchList, DNALength curPos, DNALength curGenomePos, DNALength nextGenomePos, DNALength clp, DNALength cle) {
  int i;
	cout << curPos << " " << curGenomePos << " " << nextGenomePos << " " << clp << " " << cle << endl;
  for (i = 0; i < matchList.size(); i++) {
    cout.width(8);
    cout << matchList[i].l << " ";
  }
  cout << endl;
  for (i = 0; i < matchList.size(); i++) {
    cout.width(8);
    cout << matchList[i].q << " ";
  }  
  cout << endl;
  for (i = 0; i < matchList.size(); i++) {
    cout.width(8);
    cout << matchList[i].t << " ";
  }
  cout << endl << endl;
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
														 IntervalSearchParameters &params) {
	if (pos.size() == 0) {
		return;
	}
  //
  // Search for clusters of intervals within the pos array within
  // pos[cur...next).  The value of 'next' should be the first anchor
  // outside the possible range to cluster, or the end of the anchor list.
	VectorIndex cur = 0;
	VectorIndex nPos = pos.size();
	VectorIndex next = cur + 1;
	DNALength curBoundary = 0, nextBoundary = 0;


	if (params.overlap == true) {
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
	else { 
		curBoundary = ContigStartPos(pos[cur].t);
		nextBoundary = ContigStartPos(pos[next].t);

		//
		// Advance next until the anchor is outside the interval that
		// statrts at 'cur', and is inside the same contig that the anchor
		// at cur is in.
		//

		DNALength curIntervalLength = NumRemainingBases(pos[cur].q, intervalLength);

		AdvanceIndexToPastInterval(pos, nPos, intervalLength, ContigStartPos,
															 cur, curBoundary, next, nextBoundary);


		DNALength prevStart = cur, prevEnd = next ;
		int prevSize = next - cur;
		DNALength maxStart = cur, maxEnd = next;
		int maxSize = SumAnchors(pos, cur, next);
		int curSize = maxSize;
		bool onFirst = true;
		bool recountInterval = false;


		if (curSize > minSize) {
			start.push_back(cur);
			end.push_back(next);
		}
		while ( cur < nPos ) {
			// 
			// This interval overlaps with a possible max start
			//

			if (pos[cur].t >= pos[maxStart].t and maxEnd > 0 and pos[cur].t < pos[maxEnd].t and curBoundary == nextBoundary) {
				if (curSize > maxSize) {
					maxSize = curSize;
					maxStart = cur;
					maxEnd   = next;
				}
			}
			else {
				if (maxSize > minSize) {
					//					cout << "adding " << pos[maxStart].t << " - " << pos[maxEnd].t << " " << maxSize <<endl;
					start.push_back(maxStart);
					end.push_back(maxEnd);
				}
				maxStart = cur;
				maxEnd   = next;
				maxSize  = curSize;
			}
		
			//
			// Done scoring current interval.  At this point the range
			// pos[cur...next) has been searched for a max increasing
			// interval.  Find a new range that will possibly yield a new
			// maximum interval.  
			// There are a few cases to consider:
			//
			//
			//genome  |---+----+------------+------+-----------------------|
			//  anchors  cur  cur+1        next   next+1
			//
			// Case 1.  The range on the target pos[ cur+1 ... next].t is a
			// valid interval (it is roughly the length of the read).  In this
			// case increase cur and next by 1, and search this range.
			//
			// genome  |---+----+------------+------+-----------------------|
			//            cur  cur+1        next   next+1
			// read interval   ====================
			//
			// Case 2.  The range on the target pos[cur+1 ... next] is not a
			// valid interval, and it is much longer than the length of the
			// read.  This implies that it is impossible to increase the score
			// of the read by including both 
			//
			// genome  |---+----+--------------------------------+-----+---|
			//            cur  cur+1                             next next+1
			// read interval   ==================== 
			//
			// Advance the interval until it includes the next anchor
			//
			// genome  |---+----+------------------+-------------+-----+---|
			//            cur  cur+1             cur+n          next next+1
			// read interval                     ==================== 
			// 

			// First advance pointer in anchor list.  If this advances to the
			// end, done and no need for further logic checking (break now).
		
  
			// 
			// If the next position is not within the same contig as the current,
			// advance the current to the next since it is impossible to find
			// any more intervals in the current pos.
			//
			if (curBoundary != nextBoundary) {
				cur = next;
				curBoundary = nextBoundary;

				//
				// Start the search for the first interval in the next contig
				// just after the current position.
				//
				if (next < nPos) {
					next++;
					AdvanceIndexToPastInterval(pos, nPos, intervalLength, ContigStartPos,
																		 cur, curBoundary, next, nextBoundary);
					maxSize = 0;
					maxStart = cur;
					maxEnd = next;
				}
				maxSize = curSize = SumAnchors(pos, cur, next);
			}
			else {
				//
				// The next interval is in the same contig as the current interval.
				//

				//
				// Make sure not to couble count the current interval.
				//

				int prevCur = cur;
				while (cur < next and 
							 pos[cur].t + intervalLength < pos[next].t) {
					curSize -= pos[cur].l;
					cur++;
				}
				if (cur >= nPos)
					break;

				//
				// Advance the next to outside this interval.
				//
				//pos[cur].t + intervalLength >= pos[next].t) {
				curSize += pos[next].l;
				next++;
			}

			if (next > nPos) {
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

			curBoundary = ContigStartPos(pos[cur].t);
			nextBoundary = ContigStartPos(pos[next].t);  
    
			//
			// Previously tried to advance half.  This is being removed since
			// proper heuristics are making it not necessary to use.
			//
		}

		if (curSize > minSize) {
			//			cout << "ending  " << pos[maxStart].t << " - " << pos[maxEnd].t << " " << maxSize <<endl;
			start.push_back(maxStart);
			end.push_back(maxEnd);
		}
	}
}


template<typename T_MatchList,
	       typename T_PValueFunction, 
	       typename T_WeightFunction,
      	 typename T_SequenceBoundaryDB,
         typename T_ReferenceSequence,
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
																WeightedIntervalSet &intervalQueue, T_ReferenceSequence &reference, T_Sequence &query,
																IntervalSearchParameters &params,
																vector<BasicEndpoint<ChainedMatchPos> > *chainEndpointBuffer,
																vector<Fragment> *fragmentBuffer,
                                ClusterList &clusterList,
                                VarianceAccumulator<float> &accumPValue, 
                                VarianceAccumulator<float> &accumWeight,
                                VarianceAccumulator<float> &accumNumAnchorBases,
                                const char *titlePtr=NULL
																) {

	int nRectangles;
	WeightedIntervalSet sdpiq;
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
  // pos[cur...next).  The value of 'next' should be the first anchor
  // outside the possible range to cluster, or the end of the anchor list.

	VectorIndex next = cur + 1;
	DNALength curBoundary = 0, nextBoundary = 0;
	vector<UInt> scores, prevOpt;

	vector<DNALength> start, end;
	
	
	StoreLargestIntervals(pos, ContigStartPos, intervalLength, params.minInterval, start, end, params);
	VectorIndex i;
	VectorIndex posi;
	int maxLISSize = 0;
	for (posi = 0; posi < start.size(); posi++) {
		
		lis.clear();
		lisIndices.clear();
		cur = start[posi];
		next = end[posi];
		//		cout << start[posi] << "\t" << end[posi] << "\t" << end[posi] - start[posi] << " " << pos[cur].t << " " << pos[next-1].t << endl;
		if (next - cur == 1) {
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
			//			cout << "gc: " << cur << " " << next << " " << next - cur << endl;
			if (params.globalChainType == 0) {
				lisSize = GlobalChain<ChainedMatchPos, BasicEndpoint<ChainedMatchPos> >(pos, cur, next, 
																																								lisIndices, chainEndpointBuffer);
				//				cout << lisIndices.size() << endl;
			}
			else {
        //
        //  A different call that allows for indel penalties.
        //
				//				lisSize = RestrictedGlobalChain(&pos[cur],next - cur, 0.1, lisIndices, scores, prevOpt);
				lisSize = SDPGlobalChain(&pos[cur], next-cur, lisIndices, params.minMatch, *fragmentBuffer);
				//				cout << lisIndices.size() << endl;
			}
		}
    
    // Maybe this should become a function?
		for (i = 0; i < lisIndices.size(); i++) {	lis.push_back(pos[lisIndices[i]+cur]); }

    
		// 
		// Compute pvalue of this match.
		//
		if (lis.size() > 0) {
			lisPValue = MatchPValueFunction.ComputePValue(lis, noOvpLisNBases, noOvpLisSize);
		}
		else {
			lisPValue = 0;
		}


    /*
    if (lis.size() > 0) {
      if () {
      }
    }
    */

		if (lisSize > maxLISSize) {
			maxLISSize  = lisSize;
		}

		//
		// Insert the interval into the interval queue maintaining only the 
		// top 'nBest' intervals. 
		//

		WeightedIntervalSet::iterator lastIt = intervalQueue.begin();
		MatchWeight lisWeight = MatchWeightFunction(lis);
    VectorIndex lisEnd = lis.size() - 1;

    accumPValue.Append(lisPValue);
    accumWeight.Append(lisWeight);

		if (lisPValue < params.maxPValue and lisSize > 0 and noOvpLisNBases > params.minInterval  ) {
      WeightedInterval weightedInterval(lisWeight, noOvpLisSize, noOvpLisNBases, 
                                        lis[0].t, lis[lisEnd].t + lis[lisEnd].GetLength(), 
                                        readDir, lisPValue, 
                                        lis[0].q, lis[lisEnd].q + lis[lisEnd].GetLength(), 
                                        lis);
			intervalQueue.insert(weightedInterval);
      if (weightedInterval.isOverlapping == false) {
				assert(ContigStartPos.Length(lis[0].t) == ContigStartPos.Length(lis[lis.size()-1].t));
        clusterList.Store((float)noOvpLisNBases, lis[0].t, lis[lis.size()-1].t, noOvpLisSize);
      }
		}
	}
	return maxLISSize;
}


#endif
