#ifndef COMPARE_STRINGS_H_
#define COMPARE_STRINGS_H_

template<typename T>
class DefaultCompareStrings {
 public:

	static int Compare(T lhs, T rhs) {
		return ThreeBit[lhs] - ThreeBit[rhs];
	}

	static int Compare(T *lhs, T *rhs, int length) {
		int i;
		int res;
		i = 0;
		T *lhsptr;
		T *rhsptr;
		lhsptr = lhs;
		rhsptr = rhs;
		char *lhsend = lhs + length;
		res = 0;
		while (lhsptr != lhsend and res == 0) {
			res = ThreeBit[*lhsptr] - ThreeBit[*rhsptr];
			++lhsptr;
			++rhsptr;
		}
		return res;
	}

	static int Equal(T a, T b) {
		//
		// Compare single characters.
		//
		return a == b;
	}
	static int LessThan(T *a, int aLen, T *b, int bLen) {
		int minabLen = MIN(aLen, bLen);
		if (minabLen <= 0)
			return 0;

		int cmpRes = memcmp((void*) a, (void*) b, minabLen);
		if (cmpRes < 0) {
			return 1;
		}
		else {
			return 0;
		}
	}

	static int LessThanEqual(T *a, int aLen, T *b, int bLen) {
		int minabLen = MIN(aLen, bLen);
		if (minabLen <= 0)
			return 1;
		int cmpRes = memcmp((void*) a, (void*)b, minabLen);
		if (cmpRes <= 0) {
			return 1;
		}
		else {
			return 0;
		}
	}

	static int Equal(T* a, int aLen, T *b, int bLen) {
		int minabLen = MIN(aLen, bLen);
		if (minabLen < 0)
			return 0;
		if (minabLen == 0)
			return 1;
		
		int cmpRes = memcmp((void*) a, (void*)b, minabLen);
		if (cmpRes == 0 and aLen <= bLen) {
			return 1;
		}
		else {
			return 0;
		}
	}
};


#endif
