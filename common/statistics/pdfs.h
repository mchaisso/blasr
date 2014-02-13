#ifndef P_D_F_S
#define P_D_F_S

#include <math.h>
#include "statutils.h"


float Binomial(float p, int x, int n) {
	return pow(p,n)*pow(1-p,x-n)*Choose(x,n);
}

float Poisson(float lambda, int n) {
	if (n < FactorialTableLength) {
		float lambdaToN = pow(lambda, n);
		if (lambdaToN == HUGE_VALF)
			return 1;
		return (pow(lambda, n) * exp(-lambda)) / FactorialTable[n];
	}
	else {
		return 0;
	}
}

float Gamma(float lambda, float t, int n) {
	if (n <= FactorialTableLength) {
		return (pow(lambda, n) * pow(t,n-1) * exp(-lambda*t)) / FactorialTable[n-1];
	}
	else {
		return 0;
	}
}

float Exponential(float lambda, int t) {
	return lambda * exp(-lambda*t);
}

float Normal(float mu, float sigma, float x) {
	assert(0);
	return 0;
}


#endif
