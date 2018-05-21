#ifndef PIECEWISE_MATCH_H_
#define PIECEWISE_MATCH_H_
#include <vector>
#include "datastructures/anchoring/WeightedInterval.h"
#include "algorithms/anchoring/LISPValue.h"
#include "algorithms/anchoring/LISPValueWeightor.h"
#include "datastructures/tuplelists/TupleCountTable.h" 
#include "BasicEndpoint.h"

using namespace std;
class DirMatch {
 public:
	// tpos, qpos, rotated-q pos
	unsigned int t, q, qr;
	unsigned int s: 1;
	unsigned int v: 1;
	unsigned int l: 30;
	unsigned int len;
	unsigned int score;
	DirMatch *chainPrev;	
	DirMatch() {
		t = q = s = l = 0;
		v = 1;
		chainPrev=NULL;
		score = 0;
	}
	bool operator<(const DirMatch &rhs) const {
		if (rhs.t != t) {
			return t < rhs.t;
		}
		else {
			return q < rhs.q;
		}
	}
	DirMatch&operator=(const DirMatch &rhs) {
		t = rhs.t; q = rhs.q; s = rhs.s; v= rhs.v;l=rhs.l; len=rhs.len;
	}
	DNALength GetLength() const {
		return this->len;
	}
	
	DNALength GetX() {
		return t;
	}
	DNALength GetY() {
		return q;
	}

	void SetQ(unsigned int _q) {
		q = _q;
	}

	DirMatch* GetChainPrev() {
		return chainPrev;
	}
	void SetChainPrev(DirMatch *p) {
		chainPrev = p;
	}
	
	int GetScore() {
		return score;
	}
	int SetScore(int s) {
		score = s;
		return s;
	}

	int GetXLength() {
		return len;
	}
	int GetYLength() {
		return len;
	}

};
template<typename T_MatchPos>
int ScoreIntervalSet(WeightedIntervalSet<T_MatchPos> &intervals) {
	typename WeightedIntervalSet<T_MatchPos>::iterator intv;
	int score = 0;
	for (intv = intervals.begin(); intv != intervals.end(); ++intv) {
		score += (*intv).matches.size();
	}
	return score;
}

		
template<typename T_SequenceBoundaryDB,
	typename T_PValueFunction,
	typename T_WeightFunction,
	typename T_ReferenceSequence,
	typename T_ReadSequence>
void PiecewiseMatch(  vector<ChainedMatchPos> &forwardMatches,
											vector<ChainedMatchPos> &reverseMatches,
											// How many sets to keep track of
											VectorIndex nBest, 
											
											// End search for intervals at boundary positions
											// stored in seqBoundaries
											T_SequenceBoundaryDB & seqBoundary,
											
											// First rand intervals by their p-value
											T_PValueFunction &lisPValue,  
											
											// When ranking intervals, sum over weights determined by MatchWeightFunction
											T_WeightFunction &lisWeightFn,  

											//
											// Output.
											// The increasing interval coordinates, 
											// in order by queue weight.
											WeightedIntervalSet<ChainedMatchPos> &topIntervals,
											T_ReferenceSequence &genome,
											T_ReadSequence &read,
											IntervalSearchParameters &intervalSearchParameters) {
	//
	// First try forward orientation
	//
	WeightedIntervalSet<ChainedMatchPos> forwardIntervals(nBest);
	WeightedIntervalSet<ChainedMatchPos> reverseIntervals(nBest);

	PiecewiseMatchDir(forwardMatches, reverseMatches,
										nBest, seqBoundary, lisPValue, lisWeightFn,
										forwardIntervals,
										genome, read, intervalSearchParameters);

	PiecewiseMatchDir(reverseMatches, forwardMatches, 
										nBest, seqBoundary, lisPValue, lisWeightFn,
										reverseIntervals,
										genome, read, intervalSearchParameters);



	//
	// Score the intervals.
	//
	int forwardScore = ScoreIntervalSet<ChainedMatchPos>(forwardIntervals);
	int reverseScore = ScoreIntervalSet<ChainedMatchPos>(reverseIntervals);

	typename WeightedIntervalSet<ChainedMatchPos>::iterator intv;
	for (intv = forwardIntervals.begin(); intv != forwardIntervals.end(); ++intv) {
		topIntervals.insert((WeightedInterval<ChainedMatchPos> &)(*intv));
	}

	for (intv = reverseIntervals.begin(); intv != reverseIntervals.end(); ++intv) {
		(*intv).readIndex = 1-(*intv).readIndex;
		topIntervals.insert((WeightedInterval<ChainedMatchPos> &)(*intv));
	}
}


void FilterOffDiagonalHits(vector<ChainedMatchPos> &matches, int window, int flank) {
	


}

template<typename T_SequenceBoundaryDB,
	typename T_PValueFunction,
	typename T_WeightFunction,
	typename T_ReferenceSequence,
	typename T_ReadSequence>
void PiecewiseMatchDir(  vector<ChainedMatchPos> &forwardMatches,
											vector<ChainedMatchPos> &reverseMatches,
											// How many sets to keep track of
											VectorIndex nBest, 
											
											// End search for intervals at boundary positions
											// stored in seqBoundaries
											T_SequenceBoundaryDB & seqBoundary,
											
											// First rand intervals by their p-value
											T_PValueFunction &lisPValue,  
											
											// When ranking intervals, sum over weights determined by MatchWeightFunction
											T_WeightFunction &lisWeightFn,  

											//
											// Output.
											// The increasing interval coordinates, 
											// in order by queue weight.
											WeightedIntervalSet<ChainedMatchPos> &topIntervals,
											T_ReferenceSequence &genome,
											T_ReadSequence &read,
											IntervalSearchParameters &intervalSearchParameters) {

	bool strandSwap = false;

	//
	// Transform forward+ reverse matches into forward only, but keep track of the strand of each.
	//
	vector<DirMatch> matches;
	matches.resize(forwardMatches.size() + reverseMatches.size());
 
	if (matches.size() == 0) {
		return;
	}
	int i;
	for (i = 0; i < forwardMatches.size(); i++) {
		matches[i].t = forwardMatches[i].t;
		matches[i].q = forwardMatches[i].q;
		matches[i].s = 0;
		matches[i].l = forwardMatches[i].l;
		matches[i].len = forwardMatches[i].l;
	}
	int nf=forwardMatches.size();
	for (i = 0; i < reverseMatches.size(); i++) {
		matches[i+nf].t = reverseMatches[i].t;
		//		matches[i+nf].q = read.length - reverseMatches[i].q;
		matches[i+nf].q = reverseMatches[i].q;
		matches[i+nf].s = 1;
		matches[i+nf].l = reverseMatches[i].l;
		matches[i+nf].len = reverseMatches[i].l;
	}
	
	std::sort(matches.begin(), matches.end());

	//
	// Now rotate any consecutive runs of negative strand.
	//

	int runStart, runEnd;
	runStart = runEnd = 0;
	while (runStart < matches.size()) {
		while (runStart < matches.size() and matches[runStart].s == 0) {
			runStart++;
		}
		runEnd = runStart;
		while (runEnd < matches.size() and matches[runEnd].s == 1 and matches[runEnd].q >= matches[runStart].q) {
			runEnd++;
		}
		// Rotate the run of negative strands.
		//
		if (runStart == matches.size() ) {
			break;
		}
		int boxStart = read.length - (matches[runEnd-1].q  + matches[runEnd-1].l);
		int boxEnd   = read.length - matches[runStart].q;

		for (int r = runStart; r < runEnd; r++) {
			matches[r].qr = matches[r].q;
			matches[r].q = boxStart + (matches[r].qr  - matches[runStart].qr);
			assert(matches[r].q < 4094952662);
		}
		runStart = runEnd;
	}

	int intervalLength = read.subreadEnd - read.subreadStart;
	WeightedIntervalSet<DirMatch> forwardIntervalQueue;

	vector<BasicEndpoint<DirMatch> > chainEndpointBuffer;
  vector<Fragment> fragmentBuffer;
/*
	ofstream orig("hits.orig.txt");

	for (int m = 0; m < matches.size(); m++) {
		if (matches[m].l > 40) {
			orig << matches[m].q << "\t" << matches[m].t << "\t" << matches[m].s << "\t" << matches[m].l << endl;
		}
	}
	orig.close();
*/

	FindMaxIncreasingInterval(Forward,
														matches,
														// allow for indels to stretch out the mapping of the read.
														intervalLength,
														nBest,
														seqBoundary,
														lisPValue,
														lisWeightFn,
														forwardIntervalQueue,
														genome, read,
														intervalSearchParameters,
														&chainEndpointBuffer, 
														&fragmentBuffer, read.title);
	
	WeightedIntervalSet<DirMatch>::iterator intv;
	int index = 0;
	for (intv = forwardIntervalQueue.begin();
			 intv != forwardIntervalQueue.end();
			 intv++) {
		int syntStart= 0;
		int syntEnd = 0;
		int split = 0;

		//
		// Start out removing small stretches of reverse anchors.
		//

		/*			
		exit(0);		
		*/
		//
		// Flip back any reverse hits.
		//

		for (int m = 0; m < intv->matches.size(); m++) {
			if (intv->matches[m].s == 1) {
				intv->matches[m].q = intv->matches[m].qr;
			}
		}
		
		vector<bool> remove(intv->matches.size(), false);
		while (syntEnd < intv->matches.size()) {
			//
			// Find a strand swap, move to end of strand swap
			while (syntEnd < intv->matches.size() and intv->matches[syntEnd].s == intv->matches[syntStart].s) {
				syntEnd++;
			}
			int swapEnd = syntEnd+1;
			while (swapEnd < intv->matches.size() and intv->matches[swapEnd].s == intv->matches[syntEnd].s) {
				swapEnd++;
			}
			
			if (swapEnd - syntEnd < 5) {
				int i;
				for (i = syntEnd; i < swapEnd; i++) {
					remove[i]=true;
				}
				syntEnd = swapEnd;
			}
			syntStart = syntEnd;
		}
		int i, p;
		p = 0;
		for (i = 0; i < intv->matches.size(); i++) {
			if (remove[i] == false) {
				intv->matches[p]=intv->matches[i];
				p++;
			}
		}

		intv->matches.resize(p);
/*
		ofstream mout("hits.txt");
		int m;
		for (m = 0; m < intv->matches.size(); m++) {
			mout << intv->matches[m].q << "\t" << intv->matches[m].t << "\t" << intv->matches[m].s << endl;
		}
		mout.close();
*/
		int nNeg = 0;
		for (i = 0; i< intv->matches.size(); i++) {
			if (intv->matches[i].s == 1) {
				nNeg++;
			}
		}

		syntStart = syntEnd = 0;
		/*
		ofstream mat("hits.txt");
		for (m = 0; m < intv->matches.size(); m++) {
			mat << intv->matches[m].q << "\t" << intv->matches[m].t << "\t" << intv->matches[m].s << endl;
		}
		mat.close();
		exit(0);
		*/

		int maxGap = 100000;
			
		while (syntEnd < intv->matches.size()) {
			int score = 0;
			while (syntEnd < intv->matches.size()
						 and intv->matches[syntEnd].s == intv->matches[syntStart].s) {
				score+= intv->matches[syntEnd].GetLength();
				syntEnd++;
				if (syntEnd < intv->matches.size() and
						((intv->matches[syntEnd].t - intv->matches[syntEnd-1].t >= maxGap) or
						 (intv->matches[syntEnd].q - intv->matches[syntEnd-1].q >= maxGap))) {
					break;
				}
			}
			vector<ChainedMatchPos> matchList(syntEnd-syntStart);
			int i;
			for (i = 0; i < syntEnd-syntStart;i++) {
				matchList[i].t = intv->matches[syntStart+i].t;
				matchList[i].q = intv->matches[syntStart+i].q;
				matchList[i].l = intv->matches[syntStart+i].l;
			}

			int strand = intv->matches[syntStart].s;
				 
			WeightedInterval<ChainedMatchPos> weightedInterval(score, score, score,
																												 intv->matches[syntStart].t, intv->matches[syntEnd-1].t + intv->matches[syntEnd-1].len,
																												 strand, -score,
																												 intv->matches[syntStart].q, intv->matches[syntEnd-1].q + intv->matches[syntEnd-1].len, read.length,
																												 matchList
																												 );
			topIntervals.insert(weightedInterval);			
			syntStart = syntEnd;
			split++;
		
		}
	}	
}
	


#endif
