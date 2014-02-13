#ifndef ALIGNMENT_COMPONENTS_CLUSTSER_INFOMRATION_H_
#define ALIGNMENT_COMPONENTS_CLUSTSER_INFOMRATION_H_

class ClusterInformation {
public:
  int maxClusterSize;
  float meanAnchorBasesPerRead;
  float sdAnchorBasesPerRead;
  int score;
  float pctSimilarity;
  int readLength;
  float nStdDev ;
  int numSignificant;
  int numForward, numReverse;
};


#endif
