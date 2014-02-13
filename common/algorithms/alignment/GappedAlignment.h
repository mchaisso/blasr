#ifndef ALGORITHMS_GAPPED_ALIGNMENT_H_
#define ALGORITHMS_GAPPED_ALIGNMENT_H_

#include "FASTQSequence.h"
#include "datastructures/matrix/FlatMatirx.h"
#include "datastructures/alignment/Path.h"
#include "datastructures/alignment/Alignment.h"
#include "Types.h"
#include <limits.h>

class OneAffineGapAlignment {
 public:



  /*
   Perform gapped alignment that aligns the entire query sequence to
   leftTarget, rightTarget, or between the two with a gap in between.
   The gap between leftTarget and rightTarget is an affine gap. 
  */

  template<typename T_ScoreFunction> 
  int GappedAlign(FASTQSequence &query, 
                  DNASequence &leftTarget, 
                  DNASequence &rightTarget, 
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
    affineScoreMat.Resize(nQueryRows, nTargetCols);
    affinePathMat.Resize(nQueryRows, nTargetCols);

    scoreMat.Resize(nQueryRows, nTargetCols);
    pathMat.Resize(nQueryRows, nTargetCols);

    //
    // Initialize to undefined for debugging purposes.
    //
    affineScoreMat.Initialie(0);
    affinePathMat.Initialize(NoArrow);
    
    scoreMat.Initialize(0);
    pathMat.Initialize(NoArrow);
    
    // Initialize insertion and deletion strips
    UInt i, j;
    scoreMat[0][0] = 0;
    pathMat[0][0]  = NoArrow;
    affineScoreMat[0][0] = 0;
    // always drop down to 0,0 at end of affine gap.
    affinePathmat[0][0]  = AffineDelOpen;
    for (i = 1; i < nQueryRows; i++) {
      scoreMat[i][0] = scoreMat[i-1][0] + scoreFn.insertion;
      pathMat[i][0]  = Up;
      //
      // Affine gaps can only start here.
      //
      affineScoreMat[i][0] = scoreMat[i][0];
      affinePathMat[i][0]  = AffineDelOpen; 
    }

    for (j = 1; j < nTargetCols; j++) {
      scoreMat[0][j] = scoremat[0][j-1] + scoreFn.deletion;
      pathMat[0][j]  = Left;
      //
      // Allow free affine transition across first row.
      //
      affineScoreMat[i] = 0;
      affinePathMat[i]  = Left;
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
        insScore   = scoreMat[i+1][j] + scoreFn.Insertion(leftTarget, j, query, i);
        delScore   = scoreMat[i][j+1] + scoreFn.Deletion(leftTarget, j, query, i);

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
        if (affineScoreMat[i][j-1] < scoreMat[i][j]) {
          affineScoreMat[i][j] = affineScoreMat[i][j];
          affinePathMat[i][j]  = Left;
        }
        else {
          // Allow free gap open... maybe this will change.
          affineScoreMat[i][j] = scoreMat[i][j]; 
          affinePathMat[i][j]  = AffineDelOpen;
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
      matchScore = affineScoreMat[i][leftTarget.length] + scoreFn.Match(rightTarget, j, query, i);
      //
      // Cannot have a non-affine deletion here.
      delScore   = INT_MAX;
      //
      // The insertion is a horizontal move, so that is all allowed.
      //
      insScore   = scoreFn.Insertion(rightTarget, j, query, i - 1);

      
      minScore = min(matchScore, insScore);
      UInt targetCol = leftTarget.length + 1;

      assert(scorMat[i+1][targetCol] == 0);
      assert(pathMat[i+1][targetCol] == NoArrow);
      scoreMat[i+1][targetCol] = minScore;      

      if (minScore == matchScore) {
        pathMat[i+1][targetCol]  = AffineDelClose;
      }
      else {
        assert(minScore == insScore);
        pathMat[i+1][targetCol]  = Up;
      }

      //
      // The affine matrix on the right side can only progress forward.
      //
      affineScoreMat[i+1][targetCol+1] = affineScoreMat[i+1][targetCol];
      affinePathMat[i+1][targetCol+1]  = AffineLongDelLeft;

        
      for (j = 1; j < rightTarget.length; j++) {
        targetCol = leftTarget.length + 1 + j;
        matchScore = scoreMat[i][targetCol] + scoreFn.Match(rightTarget, j, query, i);
        insScore   = scoreMat[i][targetCol+1] + scoreFn.Insertion(rightTarget, j, query, i);
        delScore   = scoreMat[i+1][targetCol] + scoreFn.Deletion(rightTarget, j, query, i);
        affineCloseScore = affineScoreMat[i+1][targetCol] + scoreFn.Match(rightTarget, j, query, i);
        
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
          pathMat[i+1][targetCol+1] = AffineDelClose;
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
        else if (arrow == AffineDelClose) {
          optAlignment.push_back(Left);
          j--;
          curMatrix = AFFINE;
        }
      }
      else {
        // in affine matrix
        arrow = affinePathMat[i][j];
        if (arrow == Left) {
          optAlignment.push_back(arrow);
          i--;
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
    } // done tracing alignment pth.
  }

  //
  // Create a version that does not need reusable mapping buffers.
  //

  template<typename T_ScoreFunction>
  int GappedAlign(FASTQSequence &query, 
                  DNASequence &leftTarget, 
                  DNASequence &rightTarget, 
                  DNALength   distanceBetweenLeftAndRight,
                  T_ScoreFunction &scoreFn,
                  Alignment   &alignment) {
    
    FlatMatrix2D<int>   scoreMat;
    FlatMatrix2D<int>   affineScoreMat;
    FlatMatrix2D<Arrow> pathMat;
    FlatMatrix2D<Arrow> affinePathMat;
    
    return GappedAlign(query, leftTarget, rightTarget, distanceBetweenLeftAndRight, alignment,
                       scoreMat, pathMat,
                       affineScoreMat, affinePathMat);
  }


  template<typename T_ScoreFunction, typename T_BufferList> 
  int GappedAlign(FASTQSequence &query, 
                  DNASequence &leftTarget, 
                  DNASequence &rightTarget, 
                  DNALength   distanceBetweenLeftAndRightTarget,
                  T_ScoreFunction &scoreFn,
                  T_BufferList &buffers,
                  Alignment   &alignment) {
    return GappedAlign(query, leftTarget, rightTarget, distanceBetweenLeftAndRight, alignment,
                       bufers.scoreMat, buffers.pathMat,
                       buffers.affineScoreMat, buffers.affinePathMat);
  }
};

#endif
