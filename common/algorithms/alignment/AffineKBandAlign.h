#ifndef AFFINE_KBAND_ALIGN_H_
#define AFFINE_KBAND_ALIGN_H_

#include <vector>
#include "../../NucConversion.h"
#include "../../defs.h"
#include "../../datastructures/matrix/FlatMatrix.h"
#include "../../datastructures/alignment/Alignment.h"
#include "KBandAlign.h"

template<typename T_QuerySequence, typename T_TargetSequence, typename T_Alignment>
int AffineKBandAlign(T_QuerySequence &pqSeq, T_TargetSequence &ptSeq,
										 int matchMat[5][5], 
										 int hpInsOpen, int hpInsExtend, int insOpen, int insExtend,
										 int del, int k,
										 vector<int> &scoreMat,
										 vector<Arrow> & pathMat,
										 vector<int> &hpInsScoreMat,
										 vector<Arrow> &hpInsPathMat,
										 vector<int> &insScoreMat,
										 vector<Arrow> &insPathMat,
										 T_Alignment &alignment, 
										 AlignmentType alignType) {

	//
	// Make a copy of the sequences that is guaranteed to be in 3-bit format 
	// for quick access to the score array.
	//									

	int INF_SCORE = INF_INT - 1000;
	T_QuerySequence qSeq;
	T_TargetSequence tSeq;
	//	CreateThreeBitSequence(pqSeq, qSeq);
	//	CreateThreeBitSequence(ptSeq, tSeq);
	qSeq.seq = pqSeq.seq;
	qSeq.length= pqSeq.length;
	tSeq.seq = ptSeq.seq;
	tSeq.length = ptSeq.length;

	
	DNALength tLen, qLen;
	SetKBoundedLengths(tSeq.length, qSeq.length, k, tLen, qLen);

	//
	//
	// Allow for length up to diagonal + k + 1 for boundary.
	// 
	// Allow for width:
	//   diagonal (1)
	//   up to k insertions (k)
	//   up to k deletions  (k)
	//   boundary on left side of matrix (1)
	// 
	//	if (qLen + k > tLen and qLen < tLen) {
	//		k = tLen - qLen +1;
	//	}
  DNALength nCols = 2*k + 1;
	VectorIndex totalMatSize = (qLen + 1) * nCols;

	
	// 
	// For now the scoreMat and path mat maintained outside this 
	// function so that they may have different sizes from the affine
	// matrices. 
	//
	
	if (scoreMat.size() < totalMatSize) {
		scoreMat.resize(totalMatSize);
		pathMat.resize(totalMatSize);
	}

	if (hpInsScoreMat.size() < totalMatSize) {
		hpInsScoreMat.resize(totalMatSize);
		hpInsPathMat.resize(totalMatSize);
		insScoreMat.resize(totalMatSize);
		insPathMat.resize(totalMatSize);
	}
	
	// 
	// Initialze matrices
	//
	std::fill(scoreMat.begin(), scoreMat.begin() + totalMatSize, 0);
	std::fill(pathMat.begin(), pathMat.begin() + totalMatSize, NoArrow);
	std::fill(hpInsScoreMat.begin(), hpInsScoreMat.begin() + totalMatSize, 0);
	std::fill(hpInsPathMat.begin(), hpInsPathMat.begin() + totalMatSize, NoArrow);
	std::fill(insScoreMat.begin(), insScoreMat.begin() + totalMatSize, 0);
	std::fill(insPathMat.begin(), insPathMat.begin() + totalMatSize, NoArrow);

	//
	// Initialize the boundaries of the DP matrix.
	//
	int q, t;
	
	if (alignType != TargetFit) {
		insScoreMat[rc2index(0, k, nCols)] = 0;
		insPathMat[rc2index(0, k, nCols)] = AffineInsOpen;
	
		for (q = 1; q <=k && q < (int) qLen + 1; q++ ){
			insScoreMat[rc2index(q, k - q, nCols)] = q * insExtend + insOpen;
			insPathMat[rc2index(q, k-q, nCols)] = AffineInsUp;
		}
	}
	else if (alignType == TargetFit) {
		// 
		// Allow free gap penalties at the beginning of the alignment.
		//
		insScoreMat[rc2index(0, k, nCols)] = 0;
		insPathMat[rc2index(0, k, nCols)]  = AffineInsOpen;
		for (q = 1; q <= k && q < (int) qLen + 1; q++ ){
			insScoreMat[rc2index(q, k - q, nCols)] = 0;
			insPathMat[rc2index(q, k-q, nCols)]    = AffineInsUp;
		}
	}
	

	//
	// Assign score for (0,0) position in matrix -- aligning a gap to a gap
	// which should just be a finished alignment.  There is no cost for
	// gap-gap alignment.
	//
	hpInsScoreMat[rc2index(0,k,nCols)] = 0;
	hpInsPathMat[rc2index(0,k,nCols)]  = AffineHPInsOpen;
	
	for (q = 1; q <= k && q < (int) qLen + 1; q++) { 
		hpInsScoreMat[rc2index(q, k - q, nCols)] = q * hpInsExtend + hpInsOpen;
		hpInsPathMat[rc2index(q,k-q,nCols)] = AffineHPInsUp;
	}

	for (t = k+1; t < (int) nCols; t++ ) {
		hpInsScoreMat[rc2index(0, t, nCols)] = INF_SCORE;//  hpInsOpen + (t - k) * hpInsExtend;//; //INF_SCORE;
		hpInsPathMat[t] = NoArrow; //AffineHPInsOpen ; //NoArrow;
		insScoreMat[t] = INF_SCORE; //insOpen + (t - k) * insExtend; //INF_SCORE;
		insPathMat[t] = NoArrow; //AffineInsOpen; //NoArrow;
	}

	for (q = 1; q <= k && q < (int) qLen + 1; q++) {
		scoreMat[rc2index(q, k - q, nCols)] = insScoreMat[rc2index(q,k-q,nCols)];
		pathMat[rc2index(q, k - q , nCols)] = AffineInsClose;
	}
	for (t = 1; t <= (int) k; t++) {
		scoreMat[rc2index(0, t + k , nCols)] = t * del;
		pathMat[rc2index(0, t + k , nCols)] = Left;
	}


	//
	// The recurrence relation here is a slight modification of the
	// standard affine gap alignment.  Deletions are non-affine.  Insertions
	// are affine with different scores for homopolymer insertions, and 
	// an affine score for mixed insertions.
	//


	int matchScore, delScore;
	int hpInsExtendScore, hpInsOpenScore, insOpenScore, insExtendScore;
	int minHpInsScore, minInsScore;
	for (q = 1; q <= (int) qLen; q++) {
		for (t = q - k; t < (int) q + k + 1; t++) {
			if (t < 1) {
				continue;
			}
			if ((DNALength) t >  tLen) {
				break;
			}

			VectorIndex upper = rc2index(q-1, k + t - q + 1, nCols);
			VectorIndex curIndex = rc2index(q, k + t - q, nCols);
			
			if (t < q + k)
				hpInsOpenScore = scoreMat[upper] + hpInsOpen;
			else
				hpInsOpenScore = INF_SCORE;
			
			//
			// The homopolymer insertion score is defined only when the previous nucleotide
			// is the same as the current, in which case the homopolymer insertion score
			// is used.  If the current and previous nucleotide in the query are different,
			// the extension is not possible, and the best that can happen is a gap open.
			//
			if (q > 1 and qSeq[q-1] == qSeq[q-2]) {
				if (t < q + k) 
 					hpInsExtendScore = hpInsScoreMat[upper] + hpInsExtend;
				else 
					hpInsExtendScore = INF_SCORE;
			}
			else {
				hpInsExtendScore = INF_SCORE;
			}
			
			//
			// Since this is only allowing insertions, this grid has only horizontal and 
			// elevation arrows.
			//
			
			if (hpInsOpenScore < hpInsExtendScore) {
				hpInsPathMat[curIndex] = AffineHPInsOpen;
				minHpInsScore = hpInsOpenScore;
			}
			else {
				hpInsPathMat[curIndex] = AffineHPInsUp;
				minHpInsScore = hpInsExtendScore;
			}

			hpInsScoreMat[curIndex] = minHpInsScore;
			if (t < q + k) {
				insOpenScore = scoreMat[upper] + insOpen;
				insExtendScore = insScoreMat[upper] + insExtend;
			}
			else {
				insOpenScore = INF_SCORE;
				insExtendScore = INF_SCORE;
			}
			
			if (insOpenScore < insExtendScore) {
				insPathMat[curIndex] = AffineInsOpen;
				minInsScore = insOpenScore;
			}
			else {
				insPathMat[curIndex] = AffineInsUp;
				minInsScore = insExtendScore;
			}
			insScoreMat[curIndex] = minInsScore;
				

			// On left boundary of k-band. 
			// do not allow deletions of t.
			if (t == q - k) {
				delScore = INF_SCORE;
			}
			else {
				// cur row = q
				// cur col = t - q 
				// prev col therefore t - q - 1
				// and offset from diagonal is k + t - q - 1
				delScore = scoreMat[rc2index(q, k + t - q - 1, nCols)] + del;
			}

			// cur row = q
			// cur col = t - q

			// cur query index = q - 1
			// cur target index = t - 1
			// therefore match row (up) = q 
			//           match col (left, but since up shifted right) = t - q
			assert(rc2index(q - 1, k + t - q, nCols) < scoreMat.size());
			assert(t-1 >= 0);
			assert(q-1 >= 0);
			matchScore = scoreMat[rc2index(q - 1, k + t - q, nCols)] + matchMat[ThreeBit[qSeq.seq[q-1]]][ThreeBit[tSeq.seq[t-1]]];

			//
			//  Possibly on right boundary of k-band, in which
			//  case do not allow insertions from q.
		
			int minScore = MIN(matchScore, MIN(delScore, MIN(minInsScore, minHpInsScore)));
			curIndex = rc2index(q, k + t - q, nCols);
			assert(curIndex < scoreMat.size());
			scoreMat[curIndex] = minScore;
			if (minScore == matchScore) {
				pathMat[curIndex] = Diagonal;
			}
			else if (minScore == delScore) {
 				pathMat[curIndex] = Left;
			}
			else if (minScore == minInsScore) {
				pathMat[curIndex] = AffineInsClose;
			}
			else {
				pathMat[curIndex] = AffineHPInsClose;
			}
		}
	}
	/*
	cout << "tracing back from: " << q << ", " << t << endl;
	cout << "match score: " << endl;
	PrintFlatMatrix(&scoreMat[0], qLen + 1, nCols, cout);
	cout << " path: " << endl;
	PrintFlatMatrix(pathMat, qLen + 1, nCols, cout);
	cout << "hp  score: " << endl;
	PrintFlatMatrix(&hpInsScoreMat[0], qLen + 1, nCols, cout);
	cout << "hp  path: " << endl;
	PrintFlatMatrix(hpInsPathMat, qLen + 1, nCols, cout);	
	cout << "normal affine ins score: " << endl;
	PrintFlatMatrix(&insScoreMat[0], qLen + 1, nCols, cout);
	cout << "normal affine ins path: " << endl;
	PrintFlatMatrix(&insPathMat[0], qLen + 1, nCols, cout);
	*/
	vector<Arrow>  optAlignment;
	// First find the end position matrix.

	int minScoreTPos, minScore;
	int minScoreQPos;
	if (alignType == Global) {
		q = qLen ;
		t = k - ((int)qLen - (int)tLen);
	}
	else if (alignType == QueryFit) {
		q = qLen;
		minScoreTPos = max(q-k,1);
    DNALength index = rc2index(qLen, k + minScoreTPos - q, nCols);
		minScore = scoreMat[index];
		for (t = q - k; t < (int) q + k + 1; t++) {
			if (t < 1) { continue;}
			if (t > tLen) { break;}
			int index = rc2index(qLen,k + t - q,nCols);
			if (scoreMat[index] < minScore) {
				minScoreTPos = t;
				minScore = scoreMat[index ];
			}
		}
		t = k - ((int)qLen - minScoreTPos);
	}
	else if (alignType == TargetFit) {
		t = tLen;

		int qStart = max(0,min((int)qLen, (int)tLen) - max(0, k - max(((int)tLen) - ((int)qLen), 0)));
		int qEnd = min(qLen, tLen + k) + 1;

		minScoreQPos = qStart;
		int index = rc2index(minScoreQPos, k - (minScoreQPos - tLen), nCols);
		minScore = scoreMat[index];
		for (q = qStart; q < qEnd; q++) {
			// add to k since this is going up.
			index = rc2index(q, k + (q - tLen), nCols);
			if (scoreMat[index] < minScore) {
				minScoreQPos = q;
				minScore     = scoreMat[index];
			}
		}
		q = minScoreQPos;
		t = (k+((int)q-(int)tLen));
	}
	
	int optScore = scoreMat[rc2index(q, t, nCols)];
	Arrow arrow;
	MatrixLabel curMatrix = Match;

	
	while ((q > 0) or
				 (q == 0 and t > k)) {
		assert(t < 2*k+1);
		if (curMatrix == Match) {
			arrow = pathMat[rc2index(q,t, nCols)];
			if (arrow == Diagonal) {
				optAlignment.push_back(arrow);
				q--;
			}
			else if (arrow == Left) {
				optAlignment.push_back(arrow);
				t--;
			}
			//
			// The following two conditions change matrices 
			// without changing coordinates, since the gap close
			// just changes state without adding to the alignment.
			//
			else if (arrow == AffineInsClose) {
				curMatrix = AffineIns;
			}
			else if (arrow == AffineHPInsClose) {
				curMatrix = AffineHPIns;
			}
		}
		else if (curMatrix == AffineHPIns) {
			//
			// The current
			arrow = hpInsPathMat[rc2index(q,t,nCols)];
			if (arrow == AffineHPInsOpen) {
				curMatrix = Match;
			}
			else if (arrow != AffineHPInsUp) {
				cout << "ERROR! Affine homopolymer insertion path matrix MUST only have UP or OPEN arrows." << endl;
				assert(0);
			}
			optAlignment.push_back(Up);
			q--;
			t++;
		}
		else if (curMatrix == AffineIns) {
			arrow = insPathMat[rc2index(q,t,nCols)];
			if (arrow == AffineInsOpen) {
				curMatrix = Match;
			}
			else if (arrow != AffineInsUp) {
				cout << "ERROR! Affine insertion path matrix MUST only have UP or OPEN arrows."<<endl;
				assert(0);
			}
			optAlignment.push_back(Up);
			q--;
			t++;
		}
		else {
			cout << "ERROR in affine local alignment, matrix is: " << curMatrix << endl;
			assert(0);
		}
	}
	//	qSeq.Free();
	//	tSeq.Free();
	std::reverse(optAlignment.begin(), optAlignment.end());
	alignment.ArrowPathToAlignment(optAlignment);
	return optScore;
}


#endif
