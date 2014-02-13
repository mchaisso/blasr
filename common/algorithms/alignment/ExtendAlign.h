#ifndef EXTEND_ALIGN_H_
#define EXTEND_ALIGN_H_
#include <vector>
#include <algorithm>
#include "KBandAlign.h"
#include "../../NucConversion.h"
#include "../../defs.h"
#include "../../datastructures/matrix/FlatMatrix.h"
#include "../../datastructures/alignment/Alignment.h"


using namespace std;
class RCToIndex {
 public:
	int qStart, tStart;
	int middleCol;
	int band;
	int nCols;

	RCToIndex() {
		qStart = 0; tStart = 0;
    band = middleCol = nCols = 0;
	}

	int operator()(int r, int c, int &index) {
		//
		// First do some error checking on the row and column to see if it
		// is within the band.
		//
		if (r < qStart) { return 0; }
		if (c < tStart) { return 0; }
		r -= qStart;
		c -= tStart;
		if (abs(r-c) > band) { return 0; } // outside band range.
		if (c < 0) { return 0; }

		if (middleCol - (r - c) >= nCols) { return  0; }
		index = (r*nCols) + (middleCol - (r - c));
		return 1;
	}
};

class BaseIndex {
 public:
	int queryPos, refPos;
	int queryAlignLength, refAlignLength;
	int QNotAtSeqBoundary(int q) {
		return q != queryAlignLength;
	}
	int TNotAtSeqBoundary(int t) {
		return t != refAlignLength;
	}	
	int QAlignLength() {
		return queryAlignLength;
	}

	int TAlignLength() {
		return refAlignLength;
	}
};

class ForwardIndex : public BaseIndex {
 public:

	int QuerySeqPos(int q) {
		return queryPos + q;
	}

	int RefSeqPos(int t) {
		return refPos + t;
	}
	
	int GetQueryStartPos(int startQ, int endQ) {
		return queryPos + startQ + 1;
	}

	int GetRefStartPos(int startT, int endT) {
		return refPos + startT + 1;
	}
	
	void OrderArrowVector(vector<Arrow> &mat) {
		reverse(mat.begin(), mat.end());
	}
};


class ReverseIndex : public BaseIndex {
 public:

	int QuerySeqPos(int q) {
		return queryPos - q;
	}

	int RefSeqPos(int t) {
		return refPos - t;
	}

	int GetQueryStartPos(int startQ, int endQ) {
		return queryPos - (endQ-1);
	}

	int GetRefStartPos(int startT, int endT) {
		return refPos - (endT-1);
	}
	
	void OrderArrowVector(vector<Arrow> &mat) {
	}

};

template<typename T_Alignment, 
         typename T_ScoreFn, 
	       typename T_QuerySeq, 
	       typename T_RefSeq, 
	       typename T_Index>
	int ExtendAlignment(T_QuerySeq &querySeq, int queryPos, 
											T_RefSeq   &refSeq,   int refPos,
											int k,
											vector<int>   &scoreMat,
											vector<Arrow> &pathMat,
											T_Alignment   &alignment,
											T_ScoreFn     &scoreFn,
											T_Index  &index,
											int minExtendNBases=1, // Require that there
											// are more than one
											// base to align.
											int maxNDrops=2 // A drop is a row where
											// the alignment is
											// extended without
											// increasing the alignment
											// score.  maxnDrops is the
											// maximum number of times
											// that one may have before
											// terminating the alignment
											// 
											) {
		//
	// Try extending an alignment in the forward direction as long the
	// maximum score that is extended is above a threshold above the
	// initial score.  This dynamically grows the alignment matrices as
	// the alignment is extended (or the limits of the alignment
	// matrices since reusable buffers are used).   
	// 

	int nCols = 2 * k + 1 + 1;  // 2*k is for search space, +1 is for the
	                            // middle band, and the last +1 is for the
                            	// boundary conditions at the beginning of
                            	// the array.

	RCToIndex rcToIndex;
	rcToIndex.band      = k;
	rcToIndex.nCols     = nCols;
	rcToIndex.middleCol = k+2-1;

	if (index.queryAlignLength  < minExtendNBases or
			index.refAlignLength < minExtendNBases) {
		//
		// One of the sequences isn't long enough to even try to extend,
		// just bail with an empty alignment.
		//
		return 0;
	}

	//
	// Preallocate arrays to be at least k long.  The full matrix may
	// not be loaded.
	//
	int matSize = nCols * (k+1);
	if (scoreMat.size() < nCols * (k+1)) {
		scoreMat.resize(nCols * (k+1));
		pathMat.resize(nCols * (k+1));
	}

	//
	// Initialize boundary conditions.
	//

	int q, t;
	// Initialize first column for insertions.
	int firstIndex;
	fill(scoreMat.begin(), scoreMat.begin() + matSize, 0);
	fill(pathMat.begin(), pathMat.begin() + matSize, NoArrow);	
	rcToIndex(0, 0, firstIndex);
	scoreMat[firstIndex] = 0;
	pathMat[firstIndex]  = NoArrow;

	// Initialize insertion penalties.
	t = 0;
	int i;
	int pi;
	for (q = 1; q <= k and index.QNotAtSeqBoundary(q-1); q++) {
		bool res = rcToIndex(q, t, i);
		assert(res);
		res = rcToIndex(q-1, t, pi);
		int qSeqPos = index.QuerySeqPos(q-1);
		scoreMat[i] = scoreMat[pi] + scoreFn.Insertion(querySeq, qSeqPos);
		pathMat[i]  = Up;
		//		cout << "initializing insertion gap penalty for " << q << " " << refPos-1 << " "  << i << " " << scoreMat[i] << endl;
	}


	// Initialize the first row for deletions.
	q = 0;
	
	for (t = 1; t <= k and index.TNotAtSeqBoundary(t-1); t++) {
		bool res = rcToIndex(q, t, i);
		assert(res);
		int previ;
		res = rcToIndex(q,t-1,previ);
		
		int qSeqPos = index.QuerySeqPos(0);
		scoreMat[i] = scoreMat[previ] + scoreFn.Deletion(querySeq, qSeqPos);
		pathMat[i]  = Left;
		//		cout << "initializing deletion gap penalty for " << ((int)queryPos)-1 << " " << t << " " << i << " " << scoreMat[i] << endl;
	}
	/*	PrintFlatMatrix(&scoreMat[0], k , nCols, cout);
	cout << endl;
	PrintFlatMatrix(&pathMat[0],  k, nCols, cout);
	cout << endl;
	*/
	int nDrops = 0;
	int prevRowMinScore = INF_INT;
	int globalMinScore = INF_INT;
	int globalMinScoreQPos = 0;
	int globalMinScoreTPos = 0;
	
	int curIndex = -1;

	int maxAlignLength = min(index.QAlignLength(), index.TAlignLength()) + maxNDrops;

	for (q = 1; (index.QNotAtSeqBoundary(q-1) and 
							 nDrops < maxNDrops and
							 q < maxAlignLength);
			 q++ ) {

		//
		// Grow the path and score matrices by another row if this has
		// extended beyond their current capacity.
		//
		if ((q+1) * nCols > scoreMat.size()) {
			scoreMat.resize((q+1)*nCols);
			pathMat.resize((q+1)*nCols);
		}

		//
		// Now score the latest row.
		//
		int curRowMinScore = INF_INT;	
		int diagLength = q;

		int tStart = max((int) 1, ((int)diagLength) - k);
		int tEnd   = min((int) (diagLength + k +1), index.TAlignLength() + 1 );
		int qSeqPos, tSeqPos;
		for (t = tStart; t < min(tEnd, maxAlignLength); t++) {
			int insIndex, delIndex, matchIndex;

			bool hasInsIndex = false, hasDelIndex = false, hasMatchIndex = false, hasCurIndex = false;

			hasCurIndex   = rcToIndex(q, t, curIndex);
			assert(hasCurIndex);

			hasDelIndex   = rcToIndex(q, t - 1, delIndex);
			hasInsIndex   = rcToIndex(q - 1, t, insIndex);
			hasMatchIndex = rcToIndex(q-1, t-1, matchIndex);

			int insScore, delScore, matchScore;
			delScore   = INF_INT;
			insScore   = INF_INT;
			matchScore = INF_INT;
			//			cout << "ins index: " << insIndex << " del: " << delIndex << " match index " << matchIndex << endl;
			qSeqPos = index.QuerySeqPos(q-1); // The offset is to allow for the boundary buffer.
			tSeqPos = index.RefSeqPos(t-1); // ditto.
			/*			if (scoreMat[insIndex] == -1) {
				cout << "bleh" << endl;
			}
			if (scoreMat[matchIndex] == -1) {
				cout << "bleh" << endl;
			}
			if (scoreMat[delIndex] == -1) {
				cout << "bleh" << endl;
			}
			
			if (scoreFn.Insertion(refSeq, (DNALength) tSeqPos, querySeq, (DNALength) qSeqPos) == -1) {
				cout << "bleh" << endl;
			}
			if (scoreFn.Deletion(refSeq, (DNALength) tSeqPos, querySeq, (DNALength) qSeqPos) == -1) {
				cout << "ugh" << endl;
			}
			if ( scoreFn.Match(refSeq, (DNALength) tSeqPos, querySeq, (DNALength) qSeqPos) == -1) {
				cout <<" gah" << endl;
				}*/

			if (hasInsIndex) {
				insScore   = scoreMat[insIndex] + scoreFn.Insertion(refSeq, (DNALength) tSeqPos, querySeq, (DNALength) qSeqPos);
			}
			if (hasDelIndex) {
				delScore   = scoreMat[delIndex] + scoreFn.Deletion(refSeq, (DNALength) tSeqPos, querySeq, (DNALength) qSeqPos);
			}
			if (hasMatchIndex) {
				matchScore = scoreMat[matchIndex] + scoreFn.Match(refSeq, (DNALength) tSeqPos, querySeq, (DNALength) qSeqPos);
			}
			/*			cout << "ins score: " << insScore << "[" << scoreMat[insIndex] << "] del score " << delScore 
					 << " [" << scoreMat[delIndex] << "] match score " << matchScore 
					 << " [" << scoreMat[matchIndex] << "] qchar " << (int) querySeq.seq[qSeqPos] << " tchar " << (int) refSeq.seq[tSeqPos] << endl;*/
			int minScore = min(matchScore, delScore);
			minScore = min(minScore, insScore);
			scoreMat[curIndex] = minScore;
			//			cout << "extend: " << qSeqPos << " " << tSeqPos << " " << minScore << endl;
			if (minScore != INF_INT) {
				if (minScore == insScore)   { pathMat[curIndex] = Up; }
				if (minScore == delScore)   { pathMat[curIndex] = Left; }
				if (minScore == matchScore) { pathMat[curIndex] = Diagonal; }
			}
			else {
				pathMat[curIndex] = NoArrow;
			}

			assert(pathMat[curIndex] != NoArrow);
			if (minScore < curRowMinScore) {
				curRowMinScore = minScore;
			}
			int nRows = q+1;
			if (minScore < globalMinScore) {
				globalMinScore = minScore;
				globalMinScoreQPos  = q;
				globalMinScoreTPos  = t;
			}

		}

		if (curRowMinScore > prevRowMinScore) {
			nDrops++;
		}
		prevRowMinScore = curRowMinScore;
	}

	int nRows = q;

	q = globalMinScoreQPos;
	t = globalMinScoreTPos;
	vector<Arrow>  optAlignment;

	rcToIndex(q,t,i);
	//
	// When the optimal score is on a cell with NoArrow, there is no
	// good alignment.  Only try and trace an alignment out if the path
	// starts on a good alignment.
	//
	if (pathMat[i] != NoArrow) {
		while(q > 0 or t > 0) {
			int res;
			res = rcToIndex(q, t, i);
			assert(res != 0);
			Arrow arrow = pathMat[i];

			optAlignment.push_back(pathMat[i]);
			if (pathMat[i] == NoArrow) {
				assert(pathMat[i] != NoArrow);
			}
			if (arrow == Diagonal) {
				q--;
				t--;
			}
			else if (arrow == Left) {
				t--;
			}
			else if (arrow == Up) {
				q--;
			}
		}
	}

	index.OrderArrowVector(optAlignment);
	alignment.ArrowPathToAlignment(optAlignment);
	alignment.qPos = index.GetQueryStartPos(q, globalMinScoreQPos);
	alignment.tPos = index.GetRefStartPos(t, globalMinScoreTPos);
	
	return globalMinScore;
}

template<typename T_Alignment, typename T_ScoreFn, typename T_QuerySeq, typename T_RefSeq>
int ExtendAlignmentForward(T_QuerySeq &querySeq, int queryPos, 
													 T_RefSeq   &refSeq,   int refPos,
													 int k,
													 vector<int>   &scoreMat,
													 vector<Arrow> &pathMat,
													 T_Alignment   &alignment,
													 T_ScoreFn     &scoreFn,
													 int minExtendNBases=1, // Require that there
																								 // are more than one
																								 // base to align.
													 int maxNDrops=2 // A drop is a row where
																					 // the alignment is
																					 // extended without
																					 // increasing the alignment
																					 // score.  maxnDrops is the
																					 // maximum number of times
																					 // that one may have before
																					 // terminating the alignment
																					 // 
													 ) {
	
	ForwardIndex forwardIndex;
	forwardIndex.queryPos = queryPos;
	forwardIndex.refPos   = refPos;
	//
	// The alignment does not include queryPos nor refPos.
	//
	forwardIndex.queryAlignLength = querySeq.length - queryPos;
	forwardIndex.refAlignLength   = refSeq.length - refPos;
	int alignScore;
	alignScore= ExtendAlignment(querySeq, queryPos, 
															refSeq, refPos, 
															k, 
															scoreMat, pathMat, 
															alignment, scoreFn, forwardIndex, minExtendNBases, maxNDrops);
	alignment.qPos = queryPos;
	alignment.tPos = refPos;
	return alignScore;
}

template<typename T_Alignment, typename T_ScoreFn, typename T_QuerySeq, typename T_RefSeq>
int ExtendAlignmentReverse(T_QuerySeq &querySeq, int queryPos, 
													 T_RefSeq   &refSeq,   int refPos,
													 int k,
													 vector<int>   &scoreMat,
													 vector<Arrow> &pathMat,
													 T_Alignment   &alignment,
													 T_ScoreFn     &scoreFn,
													 int minExtendNBases=1, // Require that there
																								 // are more than one
																								 // base to align.
													 int maxNDrops=2 // A drop is a row where
																					 // the alignment is
																					 // extended without
																					 // increasing the alignment
																					 // score.  maxnDrops is the
																					 // maximum number of times
																					 // that one may have before
																					 // terminating the alignment
																					 // 
													 ) {
	
	ReverseIndex reverseIndex;
	reverseIndex.queryPos = queryPos-1;
	reverseIndex.refPos   = refPos-1;
	reverseIndex.queryAlignLength = queryPos;
	reverseIndex.refAlignLength   = refPos;
	int alignScore;
	alignScore = ExtendAlignment(querySeq, queryPos, 
															 refSeq, refPos, 
															 k, 
															 scoreMat, pathMat, 
															 alignment, scoreFn, reverseIndex, minExtendNBases, maxNDrops);
	
	return alignScore;
	
}

											


#endif
