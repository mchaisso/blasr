#ifndef LONGEST_INCREASING_SUBSEQUENCE_H_
#define LONGEST_INCREASING_SUBSEQUENCE_H_

#include <assert.h>
#include <vector>
#include <iostream>

using namespace std;

template<typename T>
class RawValue {
	int operator()(T &t) {
		return t;
	}
};

template<typename T, typename F_IntValue>
	int BinarySearch(T *x, vector<int> &m, int i, int lenM) {
	//
	// Binary search for the largest
	//   j ≤ L such that X[M[j]] < X[i]
	//   (or set j = 0 if no such value exists)
	//

	//
	// In English: Find the longest subsequence ending
	//   at some value less than X[i].  That
  //   way, X[i] can increase this subsequence by 1.
	//
	int low = 0, high = lenM, cur;
	cur = (high + low) / 2;
	assert(cur < m.size());
	F_IntValue IntValue;
	while (low < high) {
		if (high == 1) {
			assert(cur == 0);
			return cur;
		}
		// The highest value is between cur and high.
		if (IntValue(x[m[cur]]) < IntValue(x[i])) {
			low = cur + 1;
		}
		else {
			// x[m[cur]] is above x[i], so the highest valid value is strictly below cur
			high = cur;
		}
		cur = (low + high) / 2;
	}

	// cur is either the last spot where x[m[cur]] < x[i], or
	// the first spot where x[m[cur]] > x[i];
	// Typically x[m[cur]] must be less than x[i], except on the first
	// call to this, in which case x[m[cur]] == x[i] == the first element in the array.
	if (cur == 0) {
		return cur;
	}
	else {
		return cur - 1;
	}
}


template<typename T, typename F_IntValue >
	int LongestIncreasingSubset(T *x, int xLength, vector<int> &subsetIndices, 
															vector<int> &m, vector<int> &p) {

	//
	// m[i] is the index of the LIS of length i+1
	//
	m.resize(xLength+1);
	if (xLength == 0)
		return 0;
	m[0] = -1;
	
	p.resize(xLength);
	int i;

	int maxM;
	//  On the first iteration m[0] should be set to 0.
	int lenM = 1;
	int mi;
	F_IntValue IntValue;
	for (i = 0; i < xLength; i++) { 
		//
		// Find the index of the longest increasing subsequence ending before x[i]
		//
		maxM = BinarySearch<T,F_IntValue>(x, m, i, lenM);

		//
		// p is the list of back pointers.
		// Create a reference for where the previous longest increasing subsequence
		// that x[i] is part of.
		//

		p[i] = m[maxM];

		//
		// If this creates a LIS longer than any previously known one,
		// add this to maxM.
		//
		if (maxM + 1 >= lenM) {
			assert(maxM+1 < m.size());
			m[maxM+1] = i;
			if (maxM + 1>= lenM) {
				// assume that lenM is never more than 2 greater than maxM.
				lenM++;
			}
		}
		//
		// Otherwise, x[i] cannot create a longer LIS, BUT, if 
		// it is less than the current element that ends a LIS of length maxM, it bumps
		// that element out from the slot of ending the LIS with length maxM.
		//
		else if (IntValue(x[i]) < IntValue(x[m[maxM+1]])) {
			m[maxM+1] = i;
		}
	}
	
	//
	// Trace back;
	//
	int lisLen = lenM-1;
	int lisCur = m[lisLen];
 	for (i = 0 ; i < lisLen and lisCur != -1; i++) {
	   subsetIndices.push_back(lisCur);
	   lisCur = p[lisCur];
	}
	std::reverse(subsetIndices.begin(), subsetIndices.end());	

	return lenM - 1;
}

template<typename T, typename F_IntValue>
	int LongestIncreasingSubset(T*x, int& xLength, vector<int> &subsetIndices) {
	vector<int> p;
	vector<int> m;
	return LongestIncreasingSubset<T, F_IntValue>(x, xLength, subsetIndices, m, p);
}


#endif
