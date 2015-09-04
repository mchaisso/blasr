#ifndef ALGORITHMS_ONEGAP_ALIGNMENT_H_
#define ALGORITHMS_ONEGAP_ALIGNMENT_H_

#include "FASTQSequence.h"
#include "datastructures/matrix/FlatMatrix.h"
#include "datastructures/alignment/Path.h"
#include "datastructures/alignment/Alignment.h"
#include "Types.h"
#include <limits.h>

  /*
    Perform gapped alignment that aligns the entire query sequence to
    leftTarget, rightTarget, or between the two with a gap in between.
    The gap between leftTarget and rightTarget is an affine gap. 
  */

template<typename T_QuerySequence, typename T_RefSequence, typename T_ScoreFunction> 
  int OneGapAlign(T_QuerySequence &query, 
                  T_RefSequence &leftTarget, 
                  T_RefSequence &rightTarget, 
                  DNALength   distanceBetweenLeftAndRightTarget,
                  T_ScoreFunction &scoreFn,
                  Alignment   &alignment, 
                  FlatMatrix2D<int> & scoreMat, FlatMatrix2D<Arrow> & pathMat, 
                  FlatMatrix2D<int> &affineScoreMat, FlatMatrix2D<Arrow> &affinePathMat) {

    /*
      Perform alignment that spans what is effectively two pairs of
      matrices.  This is implemented as a single matrix, however paths
      may only transition through the boundary between the two through
      the affine portion of the matrices.  

      leftTarget    rightTarget
affine   |==============|
         |   x------>x  |
         |   ^       |  |
         |===|=======|==|
         |       |
regular  |===|=+||===|=+|
         |   x  ||   |  |
         |      ||   v  |
         |======||======|

             
    */
       

  UInt nQueryRows, nLeftTargetCols, nRightTargetCols;
  UInt nTargetCols;
  nQueryRows = query.length + 1;
  nLeftTargetCols  = leftTarget.length + 1;
  nRightTargetCols = rightTarget.length;
  nTargetCols = nLeftTargetCols + nRightTargetCols;
    
  //
  // Create the matrices
  // 
  affineScoreMat.Grow(nQueryRows,  nTargetCols);
  affinePathMat.Grow(nQueryRows,  nTargetCols);

  scoreMat.Grow(nQueryRows, nTargetCols);
  pathMat.Grow(nQueryRows, nTargetCols);

  //
  // Initialize to undefined for debugging purposes.
  //
  affineScoreMat.Initialize(0);
  affinePathMat.Initialize(NoArrow);
    
  scoreMat.Initialize(0);
  pathMat.Initialize(NoArrow);
    
  // Initialize insertion and deletion strips
  UInt i, j;
  scoreMat[0][0] = 0;
  pathMat[0][0]  = NoArrow;
  affineScoreMat[0][0] = 0;
  // always drop down to 0,0 at end of affine gap.
  affinePathMat[0][0]  = AffineDelOpen;
  for (i = 1; i < nQueryRows; i++) {
    scoreMat[i][0] = scoreMat[i-1][0] + scoreFn.ins;
    pathMat[i][0]  = Up;
    //
    // Affine gaps can only start here.
    //
    affineScoreMat[i][0] = scoreMat[i][0];
    affinePathMat[i][0]  = AffineDelOpen; 
  }

  for (j = 1; j < nTargetCols; j++) {
    scoreMat[0][j] = scoreMat[0][j-1] + scoreFn.del;
    pathMat[0][j]  = Left;
    //
    // Allow free affine transition across first row.
    //
    affineScoreMat[0][j] = 0;
    affinePathMat[0][j]  = Left;
  }

    
  //  Now run the alignment.
  //  i and j index over positions in the query/target sequences, or
  //  [0,len(seq)).  Since the score mat are of length len(seq) + 1,
  //  the indices of the score matrices are [1, len(seq) + 1)
  for (i = 0; i < query.length; i++) {
      
    //
    // First align a row of the left target, allowing for transition
    // up to the big affine gap. No transitions down from affine gap
    // are allowed.
    //
    for (j = 0; j < leftTarget.length; j++) {
      //
      // First, assign the non-affine score.
      //
      int matchScore, insScore, delScore;
      // Remember mapping:
      // scoreMat[i+1][j+1] = position to fill
      // scoreMat[i+1][j] ==> back one column
      // scoreMat[i][j+1] ==> up one row.
      matchScore = scoreMat[i][j]   + scoreFn.Match(leftTarget, j, query, i);
      insScore   = scoreMat[i][j+1] + scoreFn.Insertion(leftTarget, j, query, i);
      delScore   = scoreMat[i+1][j] + scoreFn.Deletion(leftTarget, j, query, i);

      int minScore = min(matchScore, min(insScore, delScore));
      scoreMat[i+1][j+1] = minScore;
        
 
      // set path.
      if (matchScore  == minScore) {
        pathMat[i+1][j+1]  = Diagonal;
      }
      else if (insScore == minScore) {
        pathMat[i+1][j+1]  = Up;
      }
      else {
        assert(delScore == minScore);
        pathMat[i+1][j+1]  = Left;
      }

      //
      // Next, assign the affine score
      //
      if (affineScoreMat[i+1][j] < scoreMat[i+1][j+1]) {
        affineScoreMat[i+1][j+1] = affineScoreMat[i+1][j];
        affinePathMat[i+1][j+1]  = Left;
      }
      else {
        // Allow free gap open... maybe this will change.
        affineScoreMat[i+1][j+1] = scoreMat[i+1][j+1]; 
        affinePathMat[i+1][j+1]  = AffineDelOpen;
      }
    }
     
    //
    // Now align the right target, allowing a jump over the divide.  
    //
    int affineCloseScore;
    j = 0; 
    //
    // A match here may only be preceded by an affine gap close.
    //
    int matchScore, delScore, insScore, minScore;
    matchScore = affineScoreMat[i][leftTarget.length] + scoreFn.Match(rightTarget, j, query, i);
    //
    // Cannot have a non-affine deletion here.
    delScore   = INT_MAX;
    //
    // The insertion is a horizontal move, so that is all allowed.
    //
    insScore   = scoreFn.Insertion(rightTarget, j, query, i - 1);

      
    minScore = min(matchScore, insScore);
    UInt targetCol = leftTarget.length;

    assert(scoreMat[i+1][targetCol+1] == 0);
    assert(pathMat[i+1][targetCol+1] == NoArrow);
    scoreMat[i+1][targetCol+1] = minScore;      

    if (minScore == matchScore) {
      pathMat[i+1][targetCol+1]  = AffineLongDelClose;
    }
    else {
      assert(minScore == insScore);
      pathMat[i+1][targetCol+1]  = Up;
    }

    //
    // The affine matrix on the right side can only progress forward.
    //
    affineScoreMat[i+1][targetCol+1] = affineScoreMat[i+1][targetCol];
    affinePathMat[i+1][targetCol+1]  = AffineLongDelLeft;

        
    for (j = 1; j < rightTarget.length; j++) {
      targetCol = leftTarget.length + j;
      matchScore = scoreMat[i][targetCol] + scoreFn.Match(rightTarget, j, query, i);
      insScore   = scoreMat[i][targetCol+1] + scoreFn.Insertion(rightTarget, j, query, i);
      delScore   = scoreMat[i+1][targetCol] + scoreFn.Deletion(rightTarget, j, query, i);
      affineCloseScore = affineScoreMat[i][targetCol] + scoreFn.Match(rightTarget, j, query, i);
        
      minScore = min(matchScore, min(insScore, min(delScore, affineCloseScore)));
        
      scoreMat[i+1][targetCol+1] = minScore;
      if (minScore == matchScore) {
        pathMat[i+1][targetCol+1]  = Diagonal;
      }
      else if (minScore == insScore) {
        pathMat[i+1][targetCol+1] = Up;
      }
      else if (minScore == delScore) {
        pathMat[i+1][targetCol+1] = Left;
      }
      else {
        assert(minScore == affineCloseScore);
        pathMat[i+1][targetCol+1] = AffineLongDelClose;
      }

      //
      // As with before, the affine matrix on the right side can
      // only progress forward.
      // 
      affineScoreMat[i+1][targetCol+1] = affineScoreMat[i+1][targetCol];
      affinePathMat[i+1][targetCol+1]  = Left;
    } // done aligning right target
  } // done aligning full query

    //
    // Now build the alignment string
    //
  i = nQueryRows - 1;
  j = nTargetCols - 1;
  vector<Arrow> optAlignment;
    
  int REGULAR = 0;
  int AFFINE  = 1;

  int curMatrix = REGULAR;
  Arrow arrow;

  int optScore = scoreMat[i][j];
  while (i > 0 or j > 0 or curMatrix == AFFINE) {
    if (curMatrix == REGULAR) {
      arrow = pathMat[i][j];
      if (arrow == Diagonal) {
        optAlignment.push_back(arrow);
        i--;
        j--;
      }
      else if (arrow == Left) {
        optAlignment.push_back(arrow);
        
        j--;
      }
      else if (arrow == Up) {
        optAlignment.push_back(arrow);
        i--;
      }
      else if (arrow == AffineLongDelClose) {
        optAlignment.push_back(Left);
        j--;
        i--;
        curMatrix = AFFINE;
      }
    }
    else {
      // in affine matrix
      arrow = affinePathMat[i][j];
      if (arrow == Left or arrow == AffineLongDelLeft) {
        optAlignment.push_back(arrow);
        j--;
      }
      else if (arrow == AffineDelOpen) {
        //
        // no change in i nor j, and this does not result in an
        // arrow.
        // Drop down to the regular alignment matrix.
        //
        curMatrix = REGULAR;
      }
    }
    assert(arrow != NoArrow);
    //
    // Check for wrap around.
    //
    assert(i != UINT_MAX);
    assert(j != UINT_MAX);
  } // done tracing alignment path.
  std::reverse(optAlignment.begin(), optAlignment.end());    
  alignment.LongGapArrowPathToAlignment(optAlignment, distanceBetweenLeftAndRightTarget);
  return optScore;
}

//
// Create a version that does not need reusable mapping buffers.
//

template<typename T_QuerySequence, typename T_RefSequence, typename T_ScoreFunction>
  int OneGapAlign(T_QuerySequence &query, 
                  T_RefSequence &leftTarget, 
                  T_RefSequence &rightTarget, 
                  DNALength   distanceBetweenLeftAndRightTarget,
                  T_ScoreFunction &scoreFn,
                  Alignment   &alignment) {
    
  FlatMatrix2D<int>   scoreMat;
  FlatMatrix2D<int>   affineScoreMat;
  FlatMatrix2D<Arrow> pathMat;
  FlatMatrix2D<Arrow> affinePathMat;
    

  return OneGapAlign(query, leftTarget, rightTarget, 
                     distanceBetweenLeftAndRightTarget,
                     scoreFn,
                     alignment,
                     scoreMat, pathMat,
                     affineScoreMat, affinePathMat);
}


template<typename T_QuerySequence, typename T_RefSequence, typename T_ScoreFunction, typename T_BufferList> 
  int OneGapAlign(T_QuerySequence &query, 
                  T_RefSequence   &leftTarget, 
                  T_RefSequence   &rightTarget, 
                  DNALength   distanceBetweenLeftAndRightTarget,
                  T_ScoreFunction &scoreFn,
                  T_BufferList &buffers,
                  Alignment   &alignment) {

  return OneGapAlign(query, leftTarget, rightTarget, distanceBetweenLeftAndRightTarget, alignment,
                     buffers.scoreMat, buffers.pathMat,
                     buffers.affineScoreMat, buffers.affinePathMat);

}

template<typename T_QuerySequence, typename T_RefSequence, typename T_ScoreFunction, typename T_BufferList> 
  int OneGapAlign(T_QuerySequence &query,
                  T_RefSequence   &reference,
                  T_ScoreFunction &scoreFunction,
                  T_BufferList &buffers,
                  Alignment &alignment) {

  T_RefSequence leftReference, rightReference;
  UInt leftReferenceLength = min(reference.length, query.length);
  leftReference.ReferenceSubstring(reference, 0, leftReferenceLength);

  UInt rightReferenceLength = min(reference.length - leftReferenceLength, query.length);
  rightReference.ReferenceSubstring(reference, reference.length - rightReferenceLength, rightReferenceLength);
  
  DNALength distanceBetweenLeftAndRight = reference.length - rightReferenceLength - leftReferenceLength;

  assert(distanceBetweenLeftAndRight >= 0);

  return OneGapAlign(query, 
                     leftReference, 
                     rightReference, 
                     distanceBetweenLeftAndRight,
                     scoreFunction, alignment);
}



#endif
