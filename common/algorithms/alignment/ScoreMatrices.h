#ifndef SCORE_MATRICES_H_
#define SCORE_MATRICES_H_

static int QVDistanceMatrix[5][5] = {
	{-1, 1, 1, 1, 1},
	{1, -1, 1, 1, 1},
	{1, 1, -1, 1, 1},
	{1, 1, 1, -1, 1},
	{1, 1, 1, 1, 1}
};

static int EditDistanceMatrix[5][5] = {
	{0, 1, 1, 1, 1},
	{1, 0, 1, 1, 1},
	{1, 1, 0, 1, 1},
	{1, 1, 1, 0, 1},
	{1, 1, 1, 1, 1}
};

static int SMRTDistanceMatrix[5][5] = {
	{-5, 6, 6, 6, 0},
	{6, -5, 6, 6, 0},
	{6, 6, -5, 6, 0},
	{6, 6, 6, -5, 0},
	{0, 0, 0, 0,  0}
};

static int SMRTLogProbMatrix[5][5] = {
  {0, 15, 15, 15, 15},
  {15, 0, 15, 15, 15},
  {15, 15, 0, 15, 15},
  {15, 15, 15, 0, 15},
  {15, 15, 15, 15, 0},
};

static int LowMutationMatrix[5][5] = {
	{0, 5, 5, 5, 5},
	{5, 0, 5, 5, 5},
	{5, 5, 0, 5, 5},
	{5, 5, 5, 0, 5},
	{5, 5, 5,  5, 5}
};

static int LocalAlignLowMutationMatrix[5][5] = {
	{-2, 5, 5, 5, 5},
	{5, -2, 5, 5, 5},
	{5, 5, -2, 5, 5},
	{5, 5, 5, -2, 5},
	{5, 5, 5,  5, 5}
};

static int CrossMatchMatrix[5][5] = {
	{-1, 2, 2, 2, 2},
	{2, -1, 2, 2, 2},
	{2, 2, -1, 2, 2},
	{2, 2, 2, -1, 2},
	{2, 2, 2,  2, 2}
};


#endif
