#ifndef ALGORITHMS_ALIGNMENT_CLUSTER_PROBABILITY_H_
#define ALGORITHMS_ALIGNMENT_CLUSTER_PROBABILITY_H_

float ComputeAnchorProbability(float pMatch, int minAnchorLength) {
  assert(minAnchorLength >= 0);
  assert(pMatch < 1.00001 and pMatch > 0);

  int i;
  float totalProbability = 0.0;
  float pMisMatch = 1 - pMatch;
  for (i = 0; i < minAnchorLength; i++) {
    totalProbability += pMatch * pMisMatch;
    pMatch *= pMatch;
  }
  return 1 - totalProbability;
}

float ComputeExpectedAnchorLength(float pMatch, int minAnchorLength, float pAnchor) {
  int i = 0;
  for (i = 0; i < minAnchorLength; i++) {
    pMatch *= pMatch;
  }
  float pMismatch = 1 - pMatch;
  float pValueEpsilon = 0.000000001;
  float expAnchorLength = 0;
  
  while(pMatch*pMismatch > pValueEpsilon) {
    expAnchorLength += i * pMatch*pMismatch/pAnchor;
    pMatch *= pMatch;
  }
  return expAnchorLength;
}

float AnchorBasesPerRead(int readLength, float pAnchor) {
  return pAnchor * readLength;
}

float AnchorBasesPerReadSigma(float expAnchorBasesPerRead) {
  // Assume exponential distribution: 
  return sqrt(expAnchorBasesPerRead);
}


#endif
