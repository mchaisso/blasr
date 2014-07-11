#ifndef AFFINE_GUIDE_ALIGNMENT_H_
#define AFFINE_GUIDE_ALIGNMENT_H_

#include <vector>
#include <math.h>
#include <limits.h>
#include <sstream>
// local
#include "AlignmentUtils.h"
#include "sdp/SDPFragment.h"
// other common includes
#include "NucConversion.h"
#include "defs.h"
#include "datastructures/matrix/FlatMatrix.h"
#include "datastructures/alignment/Alignment.h"
#include "datastructures/anchoring/MatchPos.h"
#include "datastructures/matrix/Matrix.h"
#include "tuples/DNATuple.h"
#include "tuples/TupleList.h"
#include "utils/LogUtils.h"
#include "utils/PhredUtils.h"
#include "qvs/QualityValueVector.h"
#include "qvs/QualityValue.h"
#include "DistanceMatrixScoreFunction.h"
#include "GuidedAlign.h"

using namespace std;
#define LOWEST_LOG_VALUE  -700


template<typename QSequence, typename TSequence, typename T_ScoreFn>
  int AffineGuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  Alignment &guideAlignment,
                        T_ScoreFn &scoreFn,
                        int bandSize,
                        Alignment &alignment,
                        vector<int>    &scoreMat,
                        vector<Arrow>  &pathMat,
                        vector<double> &probMat,
                        vector<double> &optPathProbMat,
                        vector<float>  &lnSubPValueVect,
                        vector<float>  &lnInsPValueVect,
                        vector<float>  &lnDelPValueVect,
                        vector<float>  &lnMatchPValueVect,
                        AlignmentType alignType=Global, 
                        bool computeProb=false) {

	Guide guide;
	AlignmentToGuide(guideAlignment, guide, bandSize);
	StoreMatrixOffsets(guide);
	int guideSize = ComputeMatrixNElem(guide);

	//
	// Make a copy of the sequences that is guaranteed to be in 3-bit format for faster alignment.
	// (mabybe eventually reuse the qseq and tseq memory)
	//
	QSequence qSeq;
	TSequence tSeq;
	qSeq.Assign(origQSeq);
	tSeq.Assign(origTSeq);
	
	int matrixNElem = ComputeMatrixNElem(guide);
	StoreMatrixOffsets(guide);

  /*
   * The following code is useful to produce images of the dp-matrix.
   * Make sure the sequenences are less than 5kb each though.

   Matrix<float> probMatrix;
   probMatrix.Resize(qSeq.length, tSeq.length);
   probMatrix.Initialize(0);
   ofstream matrixOut;
   stringstream matrixOutNameStrm;
   matrixOutNameStrm << "probMatrix_"<< runIndex << ".dat";
   matrixOut.open(matrixOutNameStrm.str().c_str());
  */

  if (computeProb) { 
    //
    // Convert phred scale to proper ln for faster manipulation later on.
    //
    QVToLogPScale(qSeq.substitutionQV, qSeq.length, lnSubPValueVect);
    QVToLogPScale(qSeq.insertionQV,    qSeq.length, lnInsPValueVect);
    QVToLogPScale(qSeq.deletionQV,     qSeq.length, lnDelPValueVect);
    if (lnMatchPValueVect.size() < qSeq.length) {
      lnMatchPValueVect.resize(qSeq.length);
    }
    //
    // Normalize probability vectors so that the probability of transition from each cell is 1.
    //
    int i;
    for (i = 0; i < qSeq.length; i++) {
      float subSum  = LogSumOfTwo(lnSubPValueVect[i], QVToLogPScale(scoreFn.substitutionPrior)); // prior on substitution rate
      float denominator = LogSumOfThree(lnDelPValueVect[i], lnInsPValueVect[i], subSum);
      lnDelPValueVect[i] = lnDelPValueVect[i] - denominator;
      lnSubPValueVect[i]   = lnSubPValueVect[i] - denominator;
      lnInsPValueVect[i]   = lnInsPValueVect[i] - denominator;
      lnMatchPValueVect[i] = subSum - denominator;
    }

  }

  
	
	// 
	// Make sure the alignments can fit in the reused buffers.
	//
  vector<int>    affineInsScoreMat, affineDelScoreMat;
  vector<Arrow>  affineInsPathMat, affineDelPathMat;
  affineInsScoreMat.resize(matrixNElem);
  fill(affineInsScoreMat.begin(), affineInsScoreMat.end(), scoreFn.affineOpen);
  affineDelScoreMat.resize(matrixNElem);
  fill(affineDelScoreMat.begin(), affineDelScoreMat.end(), scoreFn.affineOpen);
  affineInsPathMat.resize(matrixNElem);
  fill(affineInsPathMat.begin(), affineInsPathMat.end(), NoArrow);
  affineDelPathMat.resize(matrixNElem);
  fill(affineDelPathMat.begin(), affineDelPathMat.end(), NoArrow);
  

	if (scoreMat.size() < matrixNElem) {
		scoreMat.resize(matrixNElem);
		pathMat.resize(matrixNElem);
		fill(scoreMat.begin(), scoreMat.end(), 0);
		fill(pathMat.begin(), pathMat.end(), NoArrow);
	}
  if (computeProb) {
    if (probMat.size() < matrixNElem) {
      probMat.resize(matrixNElem);
      optPathProbMat.resize(matrixNElem);
    }
  }

	// 
	// Initialze matrices.  Only initialize up to matrixNElem rather
	// than matrix.size() because the matrix.size() unnecessary space
	// may be allocated.
	//

	std::fill(scoreMat.begin(), scoreMat.begin() + matrixNElem, 0);
	std::fill(pathMat.begin(), pathMat.begin() + matrixNElem, NoArrow);
  if (computeProb) {
    std::fill(probMat.begin(), probMat.begin() + matrixNElem, 0);	
    std::fill(optPathProbMat.begin(), optPathProbMat.begin() + matrixNElem, 0);
  }
	//
	// Initialize boundary conditions.
	//
	int q, t;
	int bufferIndex;
	// start alignemnt at the beginning of the guide, and align to the
	// end of the guide.
	if (guide.size() == 0) {
		qSeq.Free();
		tSeq.Free();
		return 0;
	}
	int qStart = guide[1].q;
	int tStart = guide[1].t;
	int qEnd   = guide[guide.size()-1].q+1;
	int tEnd   = guide[guide.size()-1].t+1;

	GetBufferIndexFunctor GetBufferIndex;
	GetBufferIndex.seqRowOffset = qStart;
  GetBufferIndex.guideSize    = guide.size();
	int indicesAreValid, delIndexIsValid, insIndexIsValid, matchIndexIsValid;
  bufferIndex = -1;
	indicesAreValid = GetBufferIndex(guide, qStart-1, tStart-1, bufferIndex);
	assert(indicesAreValid);
	scoreMat[bufferIndex] = 0;
	pathMat[bufferIndex]  = NoArrow;
	int matchIndex, insIndex, delIndex, curIndex;

	//
	// Initialize deletion row.
	//
  if (computeProb) {
    probMat[0] = optPathProbMat[0] = 0;
  }
	for (t = tStart; t < tStart + guide[0].tPost; t++) {
    curIndex=-1;
		indicesAreValid = GetBufferIndex(guide, qStart-1, t, curIndex);
		if (indicesAreValid == 0 ) {
			cout << "QSeq" << endl;
			((DNASequence)origQSeq).PrintSeq(cout);
			cout << "TSeq" << endl;
			((DNASequence)origTSeq).PrintSeq(cout);
			assert(0);
		}
    delIndex = -1;
		delIndexIsValid = GetBufferIndex(guide, qStart-1, t-1, delIndex);
		if (delIndexIsValid) {
			if (alignType == Global) {
				scoreMat[curIndex] = scoreMat[delIndex] + scoreFn.del;
			}
			else if (alignType == Local) {
				scoreMat[curIndex] = 0;
			}
      affineDelScoreMat[curIndex] = scoreFn.affineOpen;
      affineDelPathMat[curIndex]  = AffineDelOpen;
      affineInsScoreMat[curIndex] = scoreFn.affineOpen;
      affineInsPathMat[curIndex]  = AffineInsOpen;
			pathMat[curIndex] = Left;
      if (computeProb) {
        if (qSeq.qual.Empty() == false) {
          optPathProbMat[curIndex] = probMat[curIndex] = probMat[delIndex] + QVToLogPScale(scoreFn.globalDeletionPrior);
        }
      }
		}
	}

	// 
	// Initialize stripe along the top of the grid.
	//
	for (q = qStart ; q < qStart + bandSize and q < qEnd; q++) {
    insIndex = -1;
		insIndexIsValid = GetBufferIndex(guide, q-1, tStart-1, // diagonal from t-start
																		 insIndex);
    curIndex = -1;
		indicesAreValid = GetBufferIndex(guide, q, tStart-1, curIndex);

		if (insIndexIsValid and indicesAreValid) {
      assert(insIndex >= 0);
      assert(curIndex >= 0);
      if (alignType == Global) {
        scoreMat[curIndex] = scoreMat[insIndex] + scoreFn.ins;
			}
			else {
				scoreMat[curIndex] = 0;
			}
      affineInsScoreMat[curIndex] = scoreFn.affineOpen;
      affineInsPathMat[curIndex]  = AffineInsOpen;
      affineDelScoreMat[curIndex] = scoreFn.affineOpen;
      affineDelPathMat[curIndex] = AffineDelOpen;
			pathMat[curIndex] = Up;
		}
	}

	int matchScore, insScore, delScore, 
    affineInsOpenScore, affineInsExtScore, 
    affineDelOpenScore, affineDelExtScore;
	
	for (q = qStart; q < qEnd; q++) {
		int qi = q - qStart + 1;
		int tp = guide[qi].t;
    curIndex = matchIndex = insIndex = delIndex = -1;

    //
    // Do some work that will help define when matchIndex and insIndex
    // may be used.  Once delIndex is computed once, it is valid for
    // all t positions.
    //
    int prevRowTEnd = -1;
    if ( qi > 0) { 
      //
      // Define the boundaries of the column which may access previously
      // computed cells with a match.
      //
      prevRowTEnd = guide[qi-1].t + guide[qi-1].tPost;
    }    
    
		for (t = tp - guide[qi].tPre ; t < guide[qi].t + guide[qi].tPost +1; t++) {

			if (q < qStart + bandSize and t == tp - guide[qi].tPre - 1) {
				// On the boundary condition, don't access the 1st element;
				t++;
				continue;
			}
			// Make sure the index is not past the end of the sequence.
			if (t < -1) continue;
			if (t >= tEnd) continue;

      //
      // No cells are available to use for insertion cost
      // computation. 
      //
      if (t > prevRowTEnd) {
        insIndex = -1;
      }
      if (t > prevRowTEnd + 1) {
        matchIndex = -1;
      }

			//
			// Find the indices in the buffer.  Since the rows are of
			// different sizes, one can't just use offsets from the buffer
			// index. 
			//

			if (GetBufferIndex(guide, q-1,t-1, matchIndex)) {
        assert(matchIndex >= 0);
				matchScore = scoreMat[matchIndex] + scoreFn.Match(tSeq, t, qSeq, q);
			}
			else {
				matchScore = INF_INT;
			}
			
			if (GetBufferIndex(guide, q-1, t, insIndex)) {
        assert(insIndex >= 0);
				insScore = scoreMat[insIndex] + scoreFn.Insertion(tSeq,(DNALength) t, qSeq, (DNALength)q);
        affineInsExtScore = affineInsScoreMat[insIndex] + scoreFn.affineExtend; // 0 extension 
			}
			else {
				insScore = INF_INT;
        affineInsExtScore = INF_INT;
			}
			if (GetBufferIndex(guide, q, t-1, delIndex)) {
        assert(delIndex >= 0);
				delScore = scoreMat[delIndex] + scoreFn.Deletion(tSeq, (DNALength) t, qSeq, (DNALength)q);
        affineDelExtScore = affineDelScoreMat[delIndex] + scoreFn.affineExtend;
			}
			else {
				delScore = INF_INT;
        affineDelExtScore = INF_INT;
			}
			//			cout << q << " " <<t << " " << matchScore << " " << insScore << " " << affineInsExtScore << " " << delScore << " " <<affineDelExtScore << endl;
			int minScore = MIN(matchScore, MIN(insScore, MIN(delScore, MIN(affineInsExtScore, affineDelExtScore))));
			int result   = GetBufferIndex(guide, q, t, curIndex);
			// This should only loop over valid cells.
			assert(result);
      assert(curIndex >= 0);
			scoreMat[curIndex] = minScore;
			if (minScore == INF_INT) {
				pathMat[curIndex] = NoArrow;
			}
			else {
				assert(result == 1);
				if (minScore == matchScore) {
					pathMat[curIndex] = Diagonal;
				}
				else if (minScore == delScore) {
					pathMat[curIndex] = Left;
				}
				else if (minScore == insScore) {
					pathMat[curIndex] = Up;
				}
        else if (minScore == affineInsExtScore) {
          pathMat[curIndex] = AffineInsClose;
        }
        else {
          assert (minScore == affineDelExtScore) ;
          pathMat[curIndex] = AffineDelClose;
        }
			}
      
      //
      // Set the penalty to initiate an affine gap here.
      //
      affineInsOpenScore = scoreMat[curIndex] + scoreFn.affineOpen;
      affineDelOpenScore = scoreMat[curIndex] + scoreFn.affineOpen;

      if (affineInsOpenScore == INF_INT and 
          affineInsExtScore == INF_INT) {
        cout << q << " " << t << endl;
        cout << "All infinity, bad things will happen." << endl;
        cout << "the score mat here is : " << scoreMat[curIndex] << " and path " << pathMat[curIndex] << endl;
        assert(0);
      }
      if (affineInsOpenScore < affineInsExtScore) {
        affineInsPathMat[curIndex] = AffineInsOpen;
        affineInsScoreMat[curIndex] = affineInsOpenScore;
      }
      else {
        affineInsPathMat[curIndex] = AffineInsUp;
        affineInsScoreMat[curIndex] = affineInsExtScore;
      }

      if (affineDelOpenScore < affineDelExtScore) {
        affineDelPathMat[curIndex] = AffineDelOpen;
        affineDelScoreMat[curIndex] = affineDelOpenScore;
      }
      else {
        affineDelPathMat[curIndex] = AffineDelLeft;
        affineDelScoreMat[curIndex] = affineDelExtScore;
      }
		}
	}		
	// Ok, for now just trace back from qend/tend
	q = qEnd-1;
	t = tEnd-1;
	vector<Arrow>  optAlignment;
	int bufferIndexIsValid;
  int curMatrix = Match;
	int nIns = 0;
	int nDel = 0;
	while(q >= qStart or t >= tStart) {
    bufferIndex = -1;
    //    cout << "backtrace: " << q << " "<< t << endl;
		bufferIndexIsValid = GetBufferIndex(guide, q, t, bufferIndex);
		assert(bufferIndexIsValid);
    assert(bufferIndex >= 0);
		Arrow arrow;
    //    cout << q << " "<< t << " " << curMatrix << " " << arrow << endl;
    if (curMatrix == Match) {
      arrow = pathMat[bufferIndex];
      if (arrow == NoArrow) {
        tSeq.ToAscii();
        qSeq.ToAscii();
        int gi;
        for (gi = 0; gi < guide.size(); gi++) {
          cout << guide[gi].q << " " << guide[gi].t << " " << guide[gi].tPre << " " << guide[gi].tPost << endl;
        }
					
        cout << "qseq: "<< endl;
        ((DNASequence)qSeq).PrintSeq(cout);
        cout << "tseq: "<< endl;
        ((DNASequence)tSeq).PrintSeq(cout);
        cout << "ERROR, this path has gone awry at " << q << " " << t << " !" << endl;
        exit(1);
      }

      if (arrow == Diagonal) {
        optAlignment.push_back(arrow);
        q--;
        t--;
      }
      else if (arrow == Up) {
        optAlignment.push_back(arrow);
        q--;
				nIns++;
      }
      else if (arrow == Left) {
        optAlignment.push_back(arrow);
        t--;
				nDel++;
      }
      else if (arrow == AffineInsClose) {
        optAlignment.push_back(Up);
        curMatrix = AffineIns;
        q--;
				nIns++;
      }
      else if (arrow == AffineDelClose) {
        t--;
				nDel++;
        optAlignment.push_back(Left);
        curMatrix = AffineDel;
      }
    }
    else if (curMatrix == AffineIns) {
      arrow = affineInsPathMat[bufferIndex];
      if (arrow == AffineInsOpen) {
        curMatrix = Match;
      }
      else if (arrow == AffineInsUp) {
        q--;
				nIns++;
        optAlignment.push_back(Up);
      }
      else {
        cout << "ERROR!  Reached arrow " << arrow << " at " << q << " " << t << " in affine ins path mat. That is bad." << endl;
        assert(0);
      }
    }
    else {
      assert(curMatrix == AffineDel);
      arrow = affineDelPathMat[bufferIndex];
      if (arrow == AffineDelOpen) {
        curMatrix = Match;
      }
      else if (arrow == AffineDelLeft) {
        t--;
        optAlignment.push_back(Left);
      }
      else {
        cout << "ERROR! Reached arrow " << arrow << " at " << q << " " << t << " in affine del mat. This is also bad." << endl;
        assert(0);
      }
    }
	}

	alignment.nCells = ComputeMatrixNElem(guide);
	std::reverse(optAlignment.begin(), optAlignment.end());
	alignment.qPos = qStart;
	alignment.tPos = tStart;
	alignment.ArrowPathToAlignment(optAlignment);
  //  StickPrintAlignment(alignment, qSeq, tSeq, cout);
  RemoveAlignmentPrefixGaps(alignment);
	int lastIndex = 0;
	tSeq.Free();
	qSeq.Free();
  lastIndex = -1;
	if (GetBufferIndex(guide, qEnd - 1, tEnd - 1, lastIndex)) {
		alignment.score = scoreMat[lastIndex];
		return scoreMat[lastIndex];
	}
	else {
		return 0;
	}
}


template<typename QSequence, typename TSequence, typename T_ScoreFn> //, typename T_BufferCache>
	int AffineGuidedAlign(QSequence &origQSeq, TSequence &origTSeq,
                        T_ScoreFn &scoreFn,
                        int bandSize,
                        int sdpIns, int sdpDel, float sdpIndelRate,
                        Alignment &alignment,
                        AlignmentType alignType=Global,
                        bool computeProb = false,
                        int sdpTupleSize= 8) {
	Alignment sdpAlignment;

  int alignScore = SDPAlign(origQSeq, origTSeq,
                            scoreFn, sdpTupleSize, 
                            sdpIns, sdpDel, sdpIndelRate,
														sdpAlignment); //, Local, false, false);

	int b;
	for (b = 0; b < sdpAlignment.blocks.size(); b++) {
		sdpAlignment.blocks[b].qPos += sdpAlignment.qPos;
		sdpAlignment.blocks[b].tPos += sdpAlignment.tPos;
	}
  sdpAlignment.tPos = 0;
  sdpAlignment.qPos = 0;

	return AffineGuidedAlign(origQSeq, origTSeq, sdpAlignment, scoreFn, bandSize, alignment,
                           // fill in optional parameters
                           alignType, computeProb);

}


//
// Use Case: No guide yet exists for this alignment, but using
// buffers.  Run SDP alignment first, then refine on the guide.
//

template<typename QSequence, typename TSequence, typename T_ScoreFn, typename T_BufferCache>
	int AffineGuidedAlign(QSequence &origQSeq, TSequence &origTSeq, 
                        T_ScoreFn &scoreFn,
                        int bandSize,
                        int sdpIns, int sdpDel, float sdpIndelRate,
                        T_BufferCache &buffers,
                        Alignment &alignment, 
                        AlignmentType alignType=Global,
                        bool computeProb = false,
                        int sdpTupleSize= 8) {
	
	Alignment sdpAlignment;

	int alignScore = SDPAlign(origQSeq, origTSeq,
                            scoreFn, sdpTupleSize, 
                            sdpIns, sdpDel, sdpIndelRate,
														sdpAlignment, buffers, Local, false, false);

	int b;
	for (b = 0; b < sdpAlignment.blocks.size(); b++) {
		sdpAlignment.blocks[b].qPos += sdpAlignment.qPos;
		sdpAlignment.blocks[b].tPos += sdpAlignment.tPos;
	}
  sdpAlignment.tPos = 0;
  sdpAlignment.qPos = 0;

	return AffineGuidedAlign(origQSeq, origTSeq, sdpAlignment, scoreFn, bandSize, buffers, alignment,
                           // fill in optional parameters
                           alignType, computeProb);
}

//
// Use case, guide exists, using buffers
//
template<typename QSequence, typename TSequence, typename T_ScoreFn, typename T_BufferCache>
  int AffineGuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  Alignment &guideAlignment,
                        T_ScoreFn &scoreFn,
                        int bandSize,
                        T_BufferCache &buffers,
                        Alignment &alignment, 
                        AlignmentType alignType=Global, 
                        bool computeProb=false) {
  return AffineGuidedAlign(origQSeq, origTSeq, guideAlignment, scoreFn, bandSize, alignment, 
                           buffers.scoreMat,
                           buffers.pathMat,
                           buffers.probMat,
                           buffers.optPathProbMat,
                           buffers.lnSubPValueMat,
                           buffers.lnInsPValueMat,
                           buffers.lnDelPValueMat,
                           buffers.lnMatchPValueMat, alignType, computeProb);
}

//
// Missing the use case for guide does not exist, and not using
// buffers.  This is just a very long function declaration, so it's
// not worth writing until it is needed.
//


//
// Use case, guide exists, but not using buffers.                   
//
template<typename QSequence, typename TSequence, typename T_ScoreFn>
  int AffineGuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  Alignment &guideAlignment,
                        T_ScoreFn &scoreFn,
                        int bandSize,
                        Alignment &alignment, 
                        AlignmentType alignType=Global, 
                        bool computeProb=false) {

  //  Make synonyms for members of the buffers class for easier typing.
  vector<int>    scoreMat;
  vector<Arrow>  pathMat;
  vector<double> probMat;
  vector<double> optPathProbMat;
  vector<float>  lnSubPValueVect;
  vector<float>  lnInsPValueVect;
  vector<float>  lnDelPValueVect;
  vector<float>  lnMatchPValueVect;

  return AffineGuidedAlign(origQSeq, origTSeq, guideAlignment, scoreFn, bandSize, alignment,
                           scoreMat,
                           pathMat,
                           probMat,
                           optPathProbMat,
                           lnSubPValueVect,
                           lnInsPValueVect,
                           lnDelPValueVect,
                           lnMatchPValueVect, alignType, computeProb);
}

                

#endif
