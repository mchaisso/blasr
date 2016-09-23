#ifndef GRAPH_PAPER_H_
#define GRAPH_PAPER_H_

#include "datastructures/alignment/Path.h"
#include "datastructures/matrix/FlatMatrix.h"


template<typename T_Point >
bool SetBounds(vector<T_Point> &points, DNALength &minPos, DNALength &maxPos, int axis ) {
  int i;
  DNALength maxRow = 0;
  if (points.size() == 0) {
    return false;
  }
  else {
    if (axis == 0) {
      minPos = maxPos = points[0].GetX();
    }
    else {
      minPos = maxPos = points[0].GetY();
    }
  }
  for (i = 1; i < points.size(); i++) {
    DNALength curPos;
    if (axis == 0) {
      curPos = points[i].GetX();
    }
    else {
      curPos = points[i].GetY();
    }
    if (curPos < minPos) {
      minPos = curPos;
    }
    else if (curPos > maxPos) {
      maxPos = curPos;
    }
  }
}

int GetIndex(DNALength pos, DNALength minPos, DNALength maxPos, int nBins) {
  assert(maxPos != minPos);
  float diff = pos - minPos;
  float len  = maxPos - minPos;
  float ratio = diff/len;
  return   min((DNALength)(nBins-1), 
               ((DNALength)(ratio * nBins)));
}

template<typename T_Point>
int GraphPaper(vector<T_Point> &points, 
               int nRows, int nCols,
               FlatMatrix2D<int> &bins,
               FlatMatrix2D<int> &scoreMat,
               FlatMatrix2D<Arrow> &pathMat,
               vector<bool> &onOptPath, 
			   int minScore,
			   int maxScore) {

  bins.Resize(nRows, nCols);
  bins.Fill(0);
  scoreMat.Resize(nRows+1, nCols+1);
  pathMat.Resize(nRows+1, nCols+1);
  scoreMat.Fill(0);
  pathMat.Fill(NoArrow);
  onOptPath.resize(points.size());
  fill(onOptPath.begin(), onOptPath.end(), false);

  DNALength xMin, xMax, yMin, yMax;
  SetBounds(points, xMin, xMax, 0); 
  xMax++; // make half-open interval.
  SetBounds(points, yMin, yMax, 1); 
  yMax++; // ditto.
  
  //
  // First set up the grid to optimize over.
  //
  int i;

  for (i = 0; i < points.size(); i++) {
    int rowIndex = GetIndex(points[i].GetX(), xMin, xMax, nRows);
    int colIndex = GetIndex(points[i].GetY(), yMin, yMax, nCols);
    bins[rowIndex][colIndex]+= points[i].length;
  }
  
	//
	// Compute the average bin size.
	//
	
	float averageBinSize = points.size() * points[0].length / (nRows*nCols  + 1);
	float indelPenalty = averageBinSize;
  //
  // Now optimize using DP.
  //

  // First handle boundary strips
	int r, c;
	for (r = 0; r < nRows; r++) {
	  for (c = 0; c < nCols;c++) {
		bins[r][c] = min(maxScore, bins[r][c]);
		if (bins[r][c] < minScore) {
		  bins[r][c] = 0;
		}
	  }
	}
	
	
  for (r = 1; r < nRows+1; r++) {
    scoreMat[r][0] = 0;
    pathMat[r][0]  = Up;
  }
  for (c = 1; c < nCols+1; c++) {
    scoreMat[0][c] = 0;
    pathMat[0][c]  = Left;
  }
  scoreMat[0][0] = 0;
	/*
  ofstream binsOut("bins.tab");
	*/
  for (r = 1; r < nRows + 1; r++) {
    for (c = 1; c < nCols + 1; c++) {
      int diagScore, leftScore, upScore;
			//	  binsOut << r << "\t" << c << "\t" << bins[r][c] << endl;
      diagScore = scoreMat[r-1][c-1] + bins[r-1][c-1];
	  leftScore = scoreMat[r][c-1] - indelPenalty;
	  upScore   = scoreMat[r-1][c] - indelPenalty;

      int optScore;
      Arrow optDir;
      if (diagScore >= leftScore and
          diagScore >= upScore) {
        optScore = diagScore ;
        optDir   = Diagonal;
      }
      else if (leftScore >= diagScore and 
               leftScore >= upScore) {
        optScore = leftScore;
        optDir   = Left;
      }
      else {
        optScore = upScore;
        optDir   = Up;
      }
      
      scoreMat[r][c] = optScore ;
      pathMat[r][c]   = optDir;
    }
  }
	//  binsOut.close();
/*
  ofstream matScores("matscores.tab");
  for (r = 1; r < nRows; r++) {
	int alnLength = min(r, nCols);
	int lastScore = scoreMat[r][nCols];
	matScores << r << "\t" << nCols << "\t" << lastScore << "\t" << float(lastScore)/alnLength << endl;
  }
  for (c = 1; c < nCols-1; c++) {
	int alnLength = min(c, nRows);
	int lastScore = scoreMat[nRows][c];
	matScores << nRows << "\t" << c << "\t" << lastScore << "\t" << float(lastScore)/alnLength << endl;
  }
  matScores.close();
 */ 
  r = nRows; c = nCols;
  while (r > 0 or c > 0) {
    Arrow dir = pathMat[r][c];
    pathMat[r][c] = Star;
    if (dir == Diagonal) { r--; c--; }
    else if (dir == Left) { c--; }
    else if (dir == Up) { r--; }
    if (r == 0 and c == 0) { break; }
  }

  //
  // Now mark inclusion/exclusion from matrix.
  //
  int nOpt = 0;
  for (i = 0; i < points.size(); i++) {
    int rowIndex = GetIndex(points[i].GetX(), xMin, xMax, nRows);
    int colIndex = GetIndex(points[i].GetY(), yMin, yMax, nCols);
    if (pathMat[rowIndex+1][colIndex+1] == Star) {
      onOptPath[i] = true;
    }
    else if (pathMat[rowIndex][colIndex+1] == Star) {
      onOptPath[i] = true;
    }
    else if (rowIndex + 2 < nRows and pathMat[rowIndex+2][colIndex+1] == Star) {
      onOptPath[i] = true;
    }
    else if (pathMat[rowIndex+1][colIndex] == Star) {
      onOptPath[i] = true;
    }
    else if (colIndex < nCols + 2  and pathMat[rowIndex+1][colIndex+2] == Star) {
      onOptPath[i] = true;
    }
    if (onOptPath[i]) {
      ++nOpt;
    }
  }
  return nOpt;
}

template<typename T_Point>
void RemoveOffOpt(vector<T_Point> &points, vector<bool> &optPath) {
  int i, c;
  for (i = 0, c = 0; i < points.size(); i++) {
    if (optPath[i]) {
      points[c] = points[i];
      c++;
    }
  }
  points.resize(c);
}

#endif
