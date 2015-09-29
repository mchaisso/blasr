#ifndef SW_ALIGN_H_
#define SW_ALIGN_H_
#include <algorithm>
#include <vector>

#include "AlignmentUtils.h"

#include "datastructures/matrix/FlatMatrix.h"
#include "datastructures/alignment/Path.h"
#include "datastructures/alignment/Alignment.h"

using namespace std;

//
// This is to be included by alignment.h, and not by itself.
//

template<typename T_QuerySequence, typename T_TargetSequence, typename T_Alignment, typename T_ScoreFn>
int SWAlign(T_QuerySequence &qSeq, T_TargetSequence &tSeq, 
						vector<int> &scoreMat,
						vector<Arrow> &pathMat, 
						T_Alignment &alignment,
						T_ScoreFn &scoreFn,
						AlignmentType alignType = Local,
						bool trustSequences = false,
						bool printMatrix = false
						) {
	VectorIndex nRows = qSeq.length + 1;
	VectorIndex nCols = tSeq.length + 1;
	
	VectorIndex totalMatSize = nRows * nCols;
	if (scoreMat.size() < totalMatSize) {
		scoreMat.resize(totalMatSize);
		pathMat.resize(totalMatSize);
	}

	if (nRows * nCols > 10000000) {
		cerr << "slow alignment " << nRows << " " << nCols << endl;
	}
	// 
	// Initialze matrices
	std::fill(scoreMat.begin(), scoreMat.begin() + totalMatSize, 0);
	std::fill(pathMat.begin(), pathMat.begin() + totalMatSize, NoArrow);

	//
	// Initialize boundary conditions.
	//
	int r, c;
	if (alignType == Global or alignType == ScoreGlobal
			or alignType == FrontAnchored or alignType == ScoreFrontAnchored) {
		//
		// Global alignments penalize gaps at the beginning of both
		// sequences.
		//
		for (c = 0; c < (int) tSeq.length + 1; c++ ){
			scoreMat[rc2index(0, c, tSeq.length + 1)] = scoreFn.del * c;
			pathMat[rc2index(0, c, tSeq.length + 1)]  = Left;
		}
		
		for (r = 0; r < (int) qSeq.length + 1; r++ ){ 
			scoreMat[rc2index(r,0, tSeq.length + 1)]  = scoreFn.ins * r;
			pathMat[rc2index(r, 0, tSeq.length + 1)]  = Up;
		}
	}
	else if (alignType == Local or alignType == ScoreLocal or alignType == LocalBoundaries
					 // end anchoring requires free gap penalties at the
					 // beginning of sequences.  
					 or alignType == EndAnchored or alignType == ScoreEndAnchored) {
		// 
		// Local alignments may shave off the beginning of either read.
		// No penalties at the starts of reads.
		//
		for (c = 0; c < (int) tSeq.length + 1; c++ ){
			scoreMat[rc2index(0, c, tSeq.length + 1)] = 0;
			pathMat[rc2index(0, c, tSeq.length + 1)]  = NoArrow;
		}
		for (r = 0; r < (int) qSeq.length + 1; r++ ){ 
			scoreMat[rc2index(r,0, tSeq.length + 1)]  = 0;
			pathMat[rc2index(r, 0, tSeq.length + 1)]  = NoArrow;
		}
	}
	else if (alignType == QueryFit or alignType == ScoreQueryFit) {
		//
		// Query fit allows free gaps at the beginning and end
		// of the target sequence.
		//
		for (c = 0; c < (int) tSeq.length + 1; c++ ){
			scoreMat[rc2index(0, c, tSeq.length + 1)] = 0;
			pathMat[rc2index(0, c, tSeq.length + 1)]  = Left;
		}
		for (r = 0; r < (int) qSeq.length + 1; r++ ){ 
			scoreMat[rc2index(r,0, tSeq.length + 1)]  = scoreFn.ins * r;
			pathMat[rc2index(r, 0, tSeq.length + 1)]  = Up;
		}
	}		
	else if (alignType == TargetFit or alignType == ScoreTargetFit) {
		//
		// Query fit allows free gaps at the beginning and end
		// of the target sequence.
		//
		for (c = 0; c < (int) tSeq.length + 1; c++ ){
			scoreMat[rc2index(0, c, tSeq.length + 1)] = scoreFn.del * c;
			pathMat[rc2index(0, c, tSeq.length + 1)]  = Left;
		}
		for (r = 0; r < (int) qSeq.length + 1; r++ ){ 
			scoreMat[rc2index(r,0, tSeq.length + 1)]  = 0;
			pathMat[rc2index(r, 0, tSeq.length + 1)]  = Up;
		}
	}		
	else if (alignType == Overlap        or alignType == ScoreOverlap or 
					 alignType == TSuffixQPrefix or alignType == ScoreTSuffixQPrefix) {
		//
		// Overlap alignments allow a gap at the beginning of the 
		// query, and at the end of the target.
		//
		for (r = 0; r < (int) qSeq.length + 1; r++ ){ 
			scoreMat[rc2index(r,0, tSeq.length + 1)]  = scoreFn.ins*r;
			pathMat[rc2index(r, 0, tSeq.length + 1)]  = Up;
		}
		for (c = 0; c < (int) tSeq.length + 1; c++ ){
			scoreMat[rc2index(0, c, tSeq.length + 1)] = 0;
			pathMat[rc2index(0, c, tSeq.length + 1)]  = Left;
		}
	}
	else if (alignType == TPrefixQSuffix or alignType == ScoreTPrefixQSuffix) {
		//
		// Overlap alignments allow a gap at the beginning of the 
		// query, and at the end of the target.
		//
		for (c = 0; c < (int) tSeq.length + 1; c++ ){
			scoreMat[rc2index(0, c, tSeq.length + 1)] = scoreFn.del * c;
			pathMat[rc2index(0, c, tSeq.length + 1)]  = Left;
		}
		for (r = 0; r < (int) qSeq.length + 1; r++ ){ 
			scoreMat[rc2index(r,0, tSeq.length + 1)]  = 0;
			pathMat[rc2index(r, 0, tSeq.length + 1)]  = Up;
		}
	}
	
	pathMat[0] = Diagonal;

	int match, qGap, tGap;
	
	//
	// Begin matrix pointers after the 
	int *matchScorePtr = &scoreMat[0];
	int *gapQScorePtr  = &scoreMat[1];
	int *gapTScorePtr  = &scoreMat[tSeq.length + 1 ];
	int *curScorePtr   = &scoreMat[tSeq.length + 2 ];
	Arrow *optPathPtr  = &pathMat[tSeq.length + 2];
	int minScore;

	int localMinScore = 0;
	int localMinRow = 0;
	int localMinCol = 0;
	for (r = 0; r < (int) qSeq.length; r++ ){
		for (c = 0; c < (int) tSeq.length; c++ ) {
			//
			// r+1, c+1 is the current row /col in the score and path mat.
			//

			//match = matchMat[TwoBit[qSeq.seq[r]]][TwoBit[tSeq.seq[c]]] + *matchScorePtr;
			//			qGap  = *gapQScorePtr + gap;
			//			tGap  = *gapTScorePtr + gap;
			match = scoreFn.Match(tSeq, c, qSeq, r) + scoreMat[rc2index(r,c,nCols)];
			qGap  = scoreMat[rc2index(r,c+1, nCols)] + scoreFn.Insertion(tSeq, r+1, qSeq, c);
			tGap  = scoreMat[rc2index(r+1,c, nCols)] + scoreFn.Deletion(tSeq, r, qSeq, c+1);
			minScore = MIN(match, MIN(qGap, tGap));
			if (minScore < localMinScore) {
				localMinScore = minScore;
				localMinRow = r;
				localMinCol = c;
			}

			if (minScore > 0 and 
					(alignType == Local or alignType == ScoreLocal or 
					 alignType == LocalBoundaries or 
					 alignType == EndAnchored or alignType == ScoreEndAnchored )) {
				*curScorePtr = 0;
				*optPathPtr  = NoArrow;
			}
			// This staement will get easier when the alignTypes are bitfields.
			// Not sure why this explicitly checks all conditions.
			else if (alignType == Local         or alignType == Global or
							 alignType == QueryFit      or alignType == Overlap or
							 alignType == TargetFit     or alignType == ScoreTargetFit or 
							 alignType == ScoreLocal    or alignType == ScoreGlobal or
							 alignType == ScoreQueryFit or alignType == ScoreOverlap or
							 alignType == FrontAnchored or alignType == ScoreFrontAnchored or 
							 alignType == EndAnchored   or alignType == ScoreEndAnchored or
							 alignType == LocalBoundaries or
							 alignType == TPrefixQSuffix or alignType == ScoreTPrefixQSuffix or
							 alignType == TSuffixQPrefix or alignType == ScoreTSuffixQPrefix ) {
				*curScorePtr = minScore;
				//		scoreMat[rc2index(r+1,c+1, tl)] = minScore;
				if (minScore == match) {
					*optPathPtr = Diagonal;
					//pathMat[rc2index(r+1,c+1,tl)] = Diagonal;
				}
				else if (minScore == qGap) {
					*optPathPtr = Up;
					//pathMat[rc2index(r+1,c+1, tl)] = Up;
				}
				else if (minScore == tGap) {
					*optPathPtr = Left;
					//pathMat[rc2index(r+1,c+1, tl)] = Left;
				}
			}
			++matchScorePtr;
			++gapTScorePtr;
			++gapQScorePtr;
			++curScorePtr;
			++optPathPtr;
		}
		// Done processing a row.  
		// This leaves the pointers starting at the first column in the next row
		// which is a boundary column. Advance one more.
		//
		++matchScorePtr;
		++gapTScorePtr;
		++gapQScorePtr;
		++curScorePtr;
		++optPathPtr;
	}
	//
	// Now trace back in the pairwise alignment.	
	// 

	// The location of the trace back depends on the type of alignment that is done.
  int minRow, minCol;
	if (alignType == Global or alignType == ScoreGlobal or 
			alignType == EndAnchored or alignType == ScoreEndAnchored ) {
		// start at bottom right of matrix.
		r = qSeq.length;
		c = tSeq.length;
		minRow = r;
		minCol = c;
	}
	else if (alignType == Local or alignType == ScoreLocal or
					 alignType == FrontAnchored or alignType == ScoreFrontAnchored or 
					 alignType == LocalBoundaries) {
		// start at cell that gives the highest score.
		r = localMinRow;
		c = localMinCol;
		minRow = r;
		minCol = c;
	}
	else if (alignType == QueryFit       or alignType == Overlap or 
					 alignType == ScoreQueryFit  or alignType == ScoreOverlap) {
		// Start at the point at the end of the target that gives the highest score, but has the
		// end query sequence alignment.
		
		r = nRows-1;
		int minScore = scoreMat[rc2index(nRows-1, 1, nCols)];
		minCol   = 1;
		for (c = 2; c < (int) nCols; c++ ) {
			if (scoreMat[rc2index(nRows-1, c, nCols)] < minScore) {
				minScore = scoreMat[rc2index(nRows-1, c, nCols)];
				minCol = c;
			}
		}
		c = minCol;
		minRow = nRows - 1;
	}
	else if (alignType == TargetFit or alignType == ScoreTargetFit) {
		// Start at the point at the end of the target that gives the highest score, but has the
		// end query sequence alignment.
		
		//
		// Always trace back from the end of the target.
		//
		minCol =  nCols-1;
		c = nCols-1;
		r = 0;
		int minScore = scoreMat[rc2index(1, nCols-1, nCols)];
		for (r = 2; r < (int) nRows; r++ ) {
			if (scoreMat[rc2index(r, nCols-1, nCols)] < minScore) {
				minScore = scoreMat[rc2index(r, nCols-1, nCols)];
				minRow = r;
			}
		}
		// store where to trace back from in the query.
		r = minRow;
	}
	else if (alignType == TSuffixQPrefix or alignType == ScoreTSuffixQPrefix) {
		// Start at the point at the end of the target that gives the highest score, but has the
		// end query sequence alignment.
		c = nCols - 1;
		r = 1;
		
		int minScore = scoreMat[rc2index(1, nCols-1, nCols)];
		minRow = 1;
		for (r = 2; r < (int) nRows; r++) {
			if (scoreMat[rc2index(r, nCols-1, nCols)] < minScore) {
				minScore = scoreMat[rc2index(r, nCols-1, nCols)];
				minRow = r;
			}
		}
		r = minRow;
		minCol = nCols - 1;
	}
	else if (alignType == TPrefixQSuffix or alignType == ScoreTPrefixQSuffix) {
		r = nRows-1;
		int minScore = scoreMat[rc2index(nRows-1, 1, nCols)];
		minCol   = 1;
		for (c = 2; c < (int) nCols; c++ ) {
			if (scoreMat[rc2index(nRows-1, c, nCols)] < minScore) {
				minScore = scoreMat[rc2index(nRows-1, c, nCols)];
				minCol = c;
			}
		}
		c = minCol;
		minRow = nRows - 1;
	}

	if (alignType != ScoreGlobal and
			alignType != ScoreLocal and
			alignType != ScoreQueryFit and
			alignType != ScoreOverlap and 
			alignType != ScoreTPrefixQSuffix and
			alignType != ScoreTSuffixQPrefix) {
		vector<Arrow>  optAlignment;
		Arrow arrow;
		while (((alignType == Global   or alignType == FrontAnchored ) and (r > 0 or c > 0)) or // global alignment stops at top corner
					 ((alignType == QueryFit or 
						 alignType == Overlap or 
						 alignType == TSuffixQPrefix) and r > 0) or 
					 ((alignType == TPrefixQSuffix) and c > 0) or 
					 (alignType == TargetFit and c > 0) or 
					 // local alignment stops at top corner -or- when new local alignment started.
					 ((alignType == Local    or alignType == EndAnchored or alignType == LocalBoundaries)
						and r > 0 and c > 0 and pathMat[r*nCols+c] != NoArrow)
					 ) {
			arrow = pathMat[rc2index(r, c, nCols)];

			//
			// When the alignment type is localBoundaries, it is not necessary to store
			// the actual alignment.  Only the starting positions and lengts will be stored.
			//
			if (alignType != LocalBoundaries) {
				optAlignment.push_back(arrow);
			}
			if (arrow == Diagonal) {
				r--;
				c--;
			}
			else if (arrow == Up) {
				r--;
			}
			else if (arrow == Left) {
				c--;
			}
		}
		// remove the boundary condition that is added for global alignment.
		if (alignType == LocalBoundaries and alignType != Local and alignType != EndAnchored and optAlignment.size() > 0) 
					optAlignment.pop_back();
		if (optAlignment.size() > 1) 
			std::reverse(optAlignment.begin(), optAlignment.end());
		if (optAlignment.size() > 0) 
			alignment.ArrowPathToAlignment(optAlignment);
	
		//
		// If running a local alignment, the alignment does not
		// explicityly encode the gaps at the beginning and ending of the
		// alignment.  These are stored in the qPos and tPos fields.
		//
		if (alignType == TSuffixQPrefix or alignType == TPrefixQSuffix) {
			alignment.qPos = r;
			alignment.tPos = c;
		}
		else if (alignType == Local or alignType == EndAnchored or alignType == LocalBoundaries) {
			alignment.qPos = r;
			alignment.tPos = c;
			alignment.qLength = localMinRow - alignment.qPos + 1;
			alignment.tLength = localMinCol - alignment.tPos + 1;
		}
		else if (alignType == QueryFit or alignType == TargetFit) {
			alignment.qPos = r;
			alignment.tPos = c;
		}

	}
	if (printMatrix) {
		PrintFlatMatrix( &scoreMat[0], qSeq.length + 1, tSeq.length + 1, cout);
		cout << endl;
		PrintFlatMatrix( &pathMat[0], qSeq.length + 1, tSeq.length + 1, cout);	
	}
	return scoreMat[rc2index(minRow, minCol, nCols)];
}


#endif
