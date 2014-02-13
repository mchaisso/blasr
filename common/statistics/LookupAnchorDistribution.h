#ifndef STATISTICS_LOOKUP_ANCHOR_DISTRIBUTION_H_
#define STATISTICS_LOOKUP_ANCHOR_DISTRIBUTION_H_

#include "AnchorDistributionTable.h"

int LookupAnchorDistribution(int readLength, int minMatchLength, int accuracy, 
                             float &mn, float &sdn, float &mnab, float &sdnab) {

  int kIndex, accIndex, lengthIndex;
  int returnValue = 0;

  // Major index is by accuracy
  if (accuracy < anchorReadAccuracies[0]) {
    returnValue = -2;
    accuracy    = anchorReadAccuracies[0];
  }
  else if (accuracy >= anchorReadAccuracies[1]) {
    returnValue = 2;
    accuracy    = anchorReadAccuracies[1] - anchorReadAccuracies[2];
  }

  accIndex = ( ((int)accuracy) - anchorReadAccuracies[0]) / anchorReadAccuracies[2];

  // middle index is by k 
  if (minMatchLength < anchorMinKValues[0]) {
    returnValue = -1; // signal too low
    minMatchLength = anchorMinKValues[0];
  }
  else if (minMatchLength >= anchorMinKValues[1]) {
    returnValue = 1; // signal too high
    minMatchLength = anchorMinKValues[1] - anchorMinKValues[2]; // max match length
  }

  kIndex = (minMatchLength - anchorMinKValues[0])/ anchorMinKValues[2];

  // last index is by read length
  if (readLength < anchorReadLengths[0]){ 
    returnValue = -3;
    readLength = anchorReadLengths[0];
  }
  else if (readLength >= anchorReadLengths[1]) {
    returnValue = 3;
    readLength = anchorReadLengths[1] - anchorReadLengths[2]; // max read length
  }

  lengthIndex = (readLength - anchorReadLengths[0]) / anchorReadLengths[2];


  int nLengths = (anchorReadLengths[1] - anchorReadLengths[0]) / anchorReadLengths[2];
  int nAccuracies = (anchorReadAccuracies[1] - anchorReadAccuracies[0]) / anchorReadAccuracies[2];
  int nAnchors = (anchorMinKValues[1] - anchorMinKValues[0]) / anchorMinKValues[2];
  int index = accIndex*(nLengths*nAnchors) + kIndex*nLengths + lengthIndex;
  
  mn = meanNumAnchors[index];
  sdn = sdNumAnchors[index];
  mnab = meanNumAnchorBases[index];
  sdnab = sdNumAnchorBases[index];

  return returnValue;
}

#endif
