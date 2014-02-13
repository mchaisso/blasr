#ifndef FRAGMENT_SORT_H_
#define FRAGMENT_SORT_H_


template<typename T_Fragment>
class LexicographicFragmentSort {
 public:
	int operator()(const T_Fragment &a, const T_Fragment &b) const {
		return a.LessThanXY(b);
	}
};


template<typename T_Fragment>
class LexicographicFragmentSortByY {
 public:
	int operator()(const T_Fragment &a, const T_Fragment &b) const {
		return a.LessThanYX(b);
	}
};




#endif
