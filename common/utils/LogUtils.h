#ifndef UTILS_SUM_OF_LOG_H_
#define UTILS_SUM_OF_LOG_H_


#define LOG_EPSILON   -30
#define LOG_EPSILON2  -8
#define LOG_EPSILON4  (logEpsilon/4.0)
#define LOG10 2.3025850929
#include <math.h>

double LogSumOfTwo(double value1, double value2) {
  //
  // value1 and value2 are in log space already.
  //
  double minValue, maxValue;
  minValue = value1, maxValue = value2;
  
  if (maxValue < minValue) {
    minValue = value2; maxValue = value1;
  }

  // convert to log10
  minValue *= LOG10;
  maxValue *= LOG10;
  
  double difference = minValue - maxValue;
  
  if (difference < LOG_EPSILON) {
    return maxValue / LOG10;
  }
  else if (difference < LOG_EPSILON2) {
    return (maxValue + exp(difference))/LOG10;
  }
  else {
    float expv = exp(difference);
    float log1pv = log1p(expv);
    return (maxValue + log1pv)/LOG10;
  }
}

double LogSumOfThree(double value1, double value2, double value3) {
  double minValue, maxValue, middleValue;
  if (value1 > value2 and value2 > value3) {
    maxValue = value1; middleValue = value2; minValue = value3;
  }
  else if (value1 > value3 and value3 > value2) {
    maxValue = value1; middleValue = value3; minValue = value2;
  }
  else if (value2 > value1 and value1 > value3) {
    maxValue = value2; middleValue = value1; minValue = value3;
  }
  else if (value2 > value3 and value3 > value1) {
    maxValue = value2; middleValue = value3; minValue = value1;
  }
  else if (value3 > value1 and value1 > value2) {
    maxValue = value3; middleValue = value1; minValue = value2;
  }
  else {
    maxValue = value3; middleValue = value2; minValue = value1;
  }
  return LogSumOfTwo(maxValue, LogSumOfTwo(middleValue, minValue));
}


#endif
