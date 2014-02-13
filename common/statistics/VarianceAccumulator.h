#ifndef STATUTILS_VARIANCE_ACCUMULATOR_H_
#define STATUTILS_VARIANCE_ACCUMULATOR_H_

template<typename T>
class VarianceAccumulator {
 public:
  
  T GetMean() {
    return ((1.0)*sumVal)/nSamples;
  }

  T GetVariance() {
    return (1.0*sumSqVal)/nSamples -  GetMean()*GetMean();
  }
  
  
  float GetNStdDev(T value) {
    T variance = GetVariance();
    T mean     = GetMean();
    if (variance > 0) {
      return fabs(value - mean)/(sqrt(variance));
    }
    else {
      return 0;
    }
  }
    
    
  int nSamples;
  T sumSqVal;
  T sumVal;
  T maxVal;
  T minVal;
  
  VarianceAccumulator<T>() {
    sumSqVal = 0;
    sumVal   = 0;
    nSamples = 0;
    maxVal = minVal = 0;
  }

  void Append(T v) {
    if (nSamples == 0) {
      maxVal = minVal = v;
    }
    else {
      if (maxVal < v) {
        maxVal = v;
      }
      if (minVal > v) {
        minVal = v;
      }
    }
    sumSqVal += v*v;
    sumVal   += v;
    nSamples++;
  }
};

#endif
