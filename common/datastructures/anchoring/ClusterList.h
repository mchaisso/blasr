#ifndef CLUSTER_LIST_H_
#define CLUSTER_LIST_H_

#include <vector>
using namespace std;

class ClusterList {
 public:
  vector<int> numBases;
  vector<int> numAnchors;
  vector<unsigned int> startPos;
  vector<int> indices;
  int lowerSizeLimit;
  int lowerSizeLimitNumAnchors;
  DNALength curp, cure; // cur - p, and cur - e.
  int curIndex;
  bool onContigStart;
  void Clear() {
    onContigStart = true;
    numBases.clear();
    startPos.clear();
    numAnchors.clear();
  }

  bool Store(int n, DNALength p, DNALength e, int b) {
    bool intervalIsEclipsed = true;
    bool intervalEclipses   = false;
    if (onContigStart == true) {
      curp = p;
      cure = e;
      intervalIsEclipsed = false;      
    }
    else {
      if (curp <= p and cure >= e) {
        intervalIsEclipsed = true;
      }
      else {
        if (p <= curp and e >= cure) {
          intervalEclipses = true;
        }
        intervalIsEclipsed = false;
      }
    }
    if (intervalIsEclipsed == false) {
      //
      // The current interval is unique: it does not eclipse any other
      // intervals and is not eclipsed by any.
      //
      if (n >= lowerSizeLimit) {
        
        if (intervalEclipses == false or onContigStart) {
          numBases.push_back(n);
          startPos.push_back(p);
          numAnchors.push_back(b);
          onContigStart = false;

          //
          // Record where the last interval is.  Since ther intervals
          // are sorted by position, only need to compare against the
          // current interval to make sure they are not overlapping
          // other intervals.
          //

          curp = p; cure = e;
        }
        else {
          //
          // The new interval eclipses the last one added.
          //
          if (n > numBases[numBases.size()-1]) {
            numBases[numBases.size()-1] = n;
            startPos[startPos.size()-1] = p;
            numAnchors[numAnchors.size()-1] = b;
            curp = p; cure = e;
          }
        }
      }
    }
    return !intervalIsEclipsed;
  }

  void ResetCoordinates() {
    onContigStart = true;
  }

  ClusterList() {
    lowerSizeLimit = 20;
    lowerSizeLimitNumAnchors = 1;
    curp = cure = 0;
    onContigStart = true;
    curIndex = 0;
  }

};

#endif
