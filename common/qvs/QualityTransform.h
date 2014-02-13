#ifndef QVS_QUALITY_TRANSFORM
#define QVS_QUALITY_TRANSFORM

#include <math.h>
#include "QualityValue.h"

/*
 * Base lookup table class for quality values.
 */

class QualityToProb {
 public:
	float prob[MAX_QUALITY_VALUE - MIN_QUALITY_VALUE + 1];
	float operator()(int index) {
		assert(index >= 0);
		assert(index <= MAX_QUALITY_VALUE);
		return prob[index];
	}
};

/* 
 * Create a lookup table for transforming from quality value
 * to p-value using Patrick Marks' low-end expand qv = -100*log10(p/(1-p))
 */
class LowEndExpandQualityTransform {
 public:
	void operator()(QualityToProb &qt) {
		int i;
		for (i = MIN_QUALITY_VALUE; i <= MAX_QUALITY_VALUE; i++) {
			float v = pow(10,i/-100.0);
			qt.prob[i] = 1 - v/(1+v);
		}
	}
};


#endif
