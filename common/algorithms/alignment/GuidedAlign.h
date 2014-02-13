#ifndef GUIDE_ALIGNMENT_H_
#define GUIDE_ALIGNMENT_H_

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

using namespace std;
#define LOWEST_LOG_VALUE  -700

#define MAX_BAND_SIZE 250

class GuideRow {
 public:
	int q, t;
	int tPre, tPost;
	unsigned int matrixOffset; // Where the center (q) is in the score
														 // and path matrices.
	int GetRowLength() {
		return tPost + tPre + 1;
	}
		
};

typedef vector<GuideRow> Guide;

class GetBufferIndexFunctor {
	// row + rowSeqOffset - 1 == seqRow
	// so, rowInSeq = -1, rowSeqOffset = 0 -> row = 0
 public:
	int seqRowOffset;
  int guideSize;
	int operator()(Guide &guide, int seqRow, int seqCol, int &index) {
    //
    // Whenever a previous index was found, just use the one in the array to the right
    //
    int prevIndex = -1;
    if (index != -1) {
      index++;
      return 1;
    }

		int row, col;
		row = seqRow + 1 - seqRowOffset;
		col = seqCol; // use the variable here for standardized naming.
		if (row < 0 or row > guideSize ) {
			return 0;
		}
		if (col <= guide[row].t and guide[row].t - col <= guide[row].tPre) {
			index = guide[row].matrixOffset - (guide[row].t - col);
      assert(prevIndex == -1 or prevIndex == index);
			return 1;
		}
		else if (col > guide[row].t and col - guide[row].t  <= guide[row].tPost) {
			index = guide[row].matrixOffset + (col - guide[row].t);
      assert(prevIndex == -1 or prevIndex == index);
			return 1;
		}
		return 0;
	}
};



int ComputeMatrixNElem(Guide &guide) {
	int totalSize = 0;
	int r;
	for (r = 0; r < guide.size(); r++) {
		totalSize += guide[r].GetRowLength();
    //    cout << r << " " << totalSize << endl;
    assert(guide[r].GetRowLength() >= 0);
	}
	return totalSize;
}

void StoreMatrixOffsets(Guide &guide) {
	int curMatrixSize = 0;
	int r;
	for (r = 0; r < guide.size(); r++) {
		guide[r].matrixOffset = guide[r].tPre + curMatrixSize;
		curMatrixSize += guide[r].GetRowLength();
	}
}


int AlignmentToGuide(Alignment &alignment, Guide &guide, int bandSize)  {
	guide.clear();
	if (alignment.size() == 0) {
		// no blocks to make guide, exit.
		return 0;
	}
	
	int tStart, tEnd, qStart, qEnd;
	int firstBlock   = 0;
	int lastBlock    = alignment.size() - 1;
	
	tStart = alignment.blocks[firstBlock].tPos;
	tEnd   = alignment.blocks[lastBlock].TEnd();
	qStart = alignment.blocks[firstBlock].qPos;
	qEnd   = alignment.blocks[lastBlock].QEnd();

	int qAlignLength = qEnd - qStart; 
	
	// Add one extra block for boundary conditions.
	guide.resize(qAlignLength+1);

	// Initilize the first (boundary condition) row.
	guide[0].t     = tStart - 1;
	guide[0].q     = qStart - 1;
	int drift = abs(tStart - qStart);	
	if (drift > bandSize) {
		guide[0].tPost = drift;
	}
	else {
		guide[0].tPost = bandSize;
	}
	guide[0].tPre  = 0;

	// The first row of the guide matches 
	int q = 0;
	int t = 0;
	int guideIndex = 1;
	int b;

	for (b = 0; b < alignment.blocks.size(); b++) {
		//
		// First add the match stored in block b, each block is a
		// diagonal, so that makes life easy.  
		//
		int bp;

		for (bp = 0; bp < alignment.blocks[b].length; bp++) {
			guide[guideIndex].t     = alignment.blocks[b].tPos + bp;
			guide[guideIndex].q     = alignment.blocks[b].qPos + bp;
			// 
			// This complicated logic is to determine how far back the band
			// should stretch.  The problem is that if the band stretches
			// back further than the previous row, it's possible for the
			// path matrix to go backwards into cells that should not be
			// touched.  
			//
			int tDiff = guide[guideIndex].t - guide[guideIndex-1].t;
			if (bp == 0) {
				guide[guideIndex].tPre  = guide[guideIndex].t - 
					(guide[guideIndex-1].t - guide[guideIndex-1].tPre); 
				guide[guideIndex].tPost = bandSize + abs(drift);
			}
			else {
				//
				// Within aligned blocks, align around the band size.
				//
				int fullLengthTPre = (guide[guideIndex].t - 
															(guide[guideIndex-1].t - 
															 guide[guideIndex-1].tPre));
				guide[guideIndex].tPre  = min(bandSize, fullLengthTPre);
				guide[guideIndex].tPost = min(MAX_BAND_SIZE, bandSize);
/*        if (guide[guideIndex].tPre > 500 or guide[guideIndex].tPost > 500) {
          cout << guideIndex << " " << guide[guideIndex].tPre << " " << guide[guideIndex].tPost << endl;
        }
        assert(guide[guideIndex].tPre >= 0);
        assert(guide[guideIndex].tPost >= 0);
        */
			}
			guideIndex++;
		}

		//
		// Now, widen k around regions where there is drift from the diagonal.
		//


		int diagonalLength;
		int qGap, tGap;

		if (b < alignment.blocks.size()-1) {
			qGap = alignment.blocks[b+1].qPos - alignment.blocks[b].QEnd();
			tGap = alignment.blocks[b+1].tPos - alignment.blocks[b].TEnd();
			// 
			// Drift is how far from the diagonal the next block starts at.
			//
			drift           = ComputeDrift(alignment.blocks[b], alignment.blocks[b+1]);
            //drift = min(drift, 100);
			diagonalLength  = min(qGap, tGap);
			
			int diagPos;
			int qPos, tPos;
			int qEnd, tEnd;
			
			qPos = alignment.blocks[b].QEnd();
			tPos = alignment.blocks[b].TEnd();

			qEnd = alignment.blocks[b+1].qPos;
			tEnd = alignment.blocks[b+1].tPos;

			for (diagPos = 0; diagPos < diagonalLength; diagPos++, tPos++, qPos++) {
				guide[guideIndex].t     = tPos;
				guide[guideIndex].q     = qPos;
				guide[guideIndex].tPre  = min(MAX_BAND_SIZE,(guide[guideIndex].t - 
																					 (guide[guideIndex-1].t - guide[guideIndex-1].tPre)));
				guide[guideIndex].tPost = min(MAX_BAND_SIZE, bandSize + abs(drift));
        if (guide[guideIndex].tPre > 500 or guide[guideIndex].tPost > 500) {
          cout << guideIndex << " " << guide[guideIndex].tPre << " " << guide[guideIndex].tPost << endl;
        } 
				++guideIndex;
			}
			
			//
			// If the query gap is shorter than target (there is a deletion
			// of the target),  the guide must be extended down the side of
			// the gap.  See the figure below.
			//
			//  *****
			//  **  
			//  * *
			//  *  *
			//  *   * // extend down from here.
		//  *   *
			//  *   * 
			//

			while (qPos < qEnd) {
				guide[guideIndex].t = tPos;
				guide[guideIndex].q = qPos;
				// move q down.
				qPos++;
				// keep tPos fixed, the guide is straight down here.
				guide[guideIndex].tPre = min(MAX_BAND_SIZE, guide[guideIndex].t - 
																		 (guide[guideIndex-1].t - guide[guideIndex-1].tPre)); //bandSize + abs(drift);

				guide[guideIndex].tPost =  min(MAX_BAND_SIZE, bandSize + abs(drift));
				guideIndex++;
			}
		}
	}
  //int i;
  //for (i = 0; i < guide.size(); i++) {
  //  guide[i].tPre = min(guide[i].tPre, 200);
  //  guide[i].tPost = min(guide[i].tPost, 200);
  // }
  return 1; // signal ok.
}



float QVToLogPScale(char qv) {
  return qv/-10.0;
}

void QVToLogPScale(QualityValueVector<unsigned char> &qualVect, int phredVectLength, vector<float> &lnVect) {
  if (phredVectLength > lnVect.size()) {
    lnVect.resize(phredVectLength);
  }
  int i;
  for (i = 0; i < phredVectLength; i++) {
    lnVect[i] = qualVect[i]/-10.0; //log(qualVect.ToProbability(i));
  }
}


template<typename QSequence, typename TSequence, typename T_ScoreFn>
  int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  Alignment &guideAlignment,
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
    assert(matrixNElem >= 0);
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
			pathMat[curIndex] = Up;
      if (computeProb) {
        if (qSeq.qual.Empty() == false) {
          optPathProbMat[curIndex] = probMat[curIndex] = probMat[insIndex] + QVToLogPScale(scoreFn.Insertion(tSeq,(DNALength) 0, qSeq, (DNALength)q));
        }
      }
		}
	}

	int matchScore, insScore, delScore;
	
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
			}
			else {
				insScore = INF_INT;
			}
			if (GetBufferIndex(guide, q, t-1, delIndex)) {
        assert(delIndex >= 0);
				delScore = scoreMat[delIndex] + scoreFn.Deletion(tSeq, (DNALength) t, qSeq, (DNALength)q);
			}
			else {
				delScore = INF_INT;
			}
			
			int minScore = MIN(matchScore, MIN(insScore, delScore));
			int result   = GetBufferIndex(guide, q, t, curIndex);
			// This should only loop over valid cells.
			assert(result);
      assert(curIndex >= 0);
			scoreMat[curIndex] = minScore;
			if (minScore == INF_INT) {
				pathMat[curIndex] = NoArrow;
        if (computeProb) {
          probMat[curIndex] = 1;
          optPathProbMat[curIndex] = 0;
        }
			}
			else {
				assert(result == 1);
				if (minScore == matchScore) {
					pathMat[curIndex] = Diagonal;
				}
				else if (minScore == delScore) {
					pathMat[curIndex] = Left;
				}
				else {
					pathMat[curIndex] = Up;
				}

				float pMisMatch, pIns, pDel;
				// Assign these to anything over 1 to signal they are not assigned.
				pMisMatch = 2;
				pIns      = 2;
				pDel      = 2;
        if (computeProb) {
          if (matchScore != INF_INT) {
            pMisMatch = QVToLogPScale(scoreFn.NormalizedMatch(tSeq, t, qSeq, q));
          }
          if (insScore != INF_INT) {
            pIns = QVToLogPScale(scoreFn.NormalizedInsertion(tSeq, t, qSeq, q));
          }
          if (delScore != INF_INT) {
            pDel = QVToLogPScale(scoreFn.NormalizedDeletion(tSeq, t, qSeq, q));
          }

          if (qSeq.qual.Empty() == false) {
            if (matchScore != INF_INT and delScore != INF_INT and insScore != INF_INT) {
              probMat[curIndex] = LogSumOfThree(probMat[matchIndex] + pMisMatch,
                                                probMat[delIndex] + pDel,
                                                probMat[insIndex] + pIns);
            }
            else if (matchScore != INF_INT and delScore != INF_INT) {
              probMat[curIndex] = LogSumOfTwo(probMat[matchIndex] + pMisMatch,
                                              probMat[delIndex] + pDel);
            }
            else if (matchScore != INF_INT and insScore != INF_INT) {
              probMat[curIndex] = LogSumOfTwo(probMat[matchIndex] + pMisMatch,
                                              probMat[insIndex] + pIns);					
            }
            else if (insScore != INF_INT and delScore != INF_INT) {
              probMat[curIndex] = LogSumOfTwo(probMat[delIndex] + pDel,
                                              probMat[insIndex] + pIns);
            }
            else if (matchScore != INF_INT) {
              probMat[curIndex] = probMat[matchIndex] + pMisMatch;
            }
            else if (delScore != INF_INT) {
              probMat[curIndex] = probMat[delIndex] + pDel;
            }
            else if (insScore != INF_INT) {
              probMat[curIndex] = probMat[insIndex] + pIns;
            }
            //
            // Not normalizing probabilities, but using value as if it
            // was a probability later on, so cap at 0 (= log 1).
            //
            if (probMat[curIndex] > 0) {
              probMat[curIndex] = 0;
            }
            assert(!isnan(probMat[curIndex]));
					}
				}
			}
		}
	}		
	// Ok, for now just trace back from qend/tend
	q = qEnd-1;
	t = tEnd-1;
	vector<Arrow>  optAlignment;
	int bufferIndexIsValid;
	while(q >= qStart or t >= tStart) {
    bufferIndex = -1;
		bufferIndexIsValid = GetBufferIndex(guide, q, t, bufferIndex);
		assert(bufferIndexIsValid);
    assert(bufferIndex >= 0);
		Arrow arrow;
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
		optAlignment.push_back(arrow);
		if (arrow == Diagonal) {
			q--;
			t--;
		}
		else if (arrow == Up) {
			q--;
		}
		else if (arrow == Left) {
			t--;
		}
	}

	alignment.nCells = ComputeMatrixNElem(guide);
	std::reverse(optAlignment.begin(), optAlignment.end());
	alignment.qPos = qStart;
	alignment.tPos = tStart;
	alignment.ArrowPathToAlignment(optAlignment);
  RemoveAlignmentPrefixGaps(alignment);
	int lastIndex = 0;
	tSeq.Free();
	qSeq.Free();
  lastIndex = -1;
	if (GetBufferIndex(guide, qEnd - 1, tEnd - 1, lastIndex)) {
		alignment.score = scoreMat[lastIndex];
    if (computeProb) {
      alignment.probScore = probMat[lastIndex];
    }
		return scoreMat[lastIndex];
	}
	else {
		return 0;
	}
}


template<typename QSequence, typename TSequence, typename T_ScoreFn> //, typename T_BufferCache>
	int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq,
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

	return GuidedAlign(origQSeq, origTSeq, sdpAlignment, scoreFn, bandSize, alignment,
                     // fill in optional parameters
                     alignType, computeProb);

}


//
// Use Case: No guide yet exists for this alignment, but using
// buffers.  Run SDP alignment first, then refine on the guide.
//

template<typename QSequence, typename TSequence, typename T_ScoreFn, typename T_BufferCache>
	int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq, 
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

	return GuidedAlign(origQSeq, origTSeq, sdpAlignment, scoreFn, bandSize, buffers, alignment,
                     // fill in optional parameters
                     alignType, computeProb);
}

//
// Use case, guide exists, using buffers
//
template<typename QSequence, typename TSequence, typename T_ScoreFn, typename T_BufferCache>
int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  Alignment &guideAlignment,
								T_ScoreFn &scoreFn,
                int bandSize,
                T_BufferCache &buffers,
								Alignment &alignment, 
								AlignmentType alignType=Global, 
                bool computeProb=false) {
  return GuidedAlign(origQSeq, origTSeq, guideAlignment, scoreFn, bandSize, alignment, 
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
int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  Alignment &guideAlignment,
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

  return GuidedAlign(origQSeq, origTSeq, guideAlignment, scoreFn, bandSize, alignment,
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
