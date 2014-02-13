#ifndef SIMULATOR_CDF_MAP_H_
#define SIMULATOR_CDF_MAP_H_
#include "statistics/statutils.h"
#include <algorithm>

template<typename T_Data>
class CDFMap {
 public:
	vector<int> cdf;
	vector<T_Data> data;
    /* Tests:
     * indices                0 1 2 3 4  5  6
     * lengths:               1 3 5 9 10 10 11
     * lengthHistogram.data:  1 3 5 9 10 11
     * lengthHistogram.cdf :  1 2 3 4 6  7
     *
     * lengths:               1 3 5 9 10 11
     * lengthHistogram.data:  1 3 5 9 10 11
     * lengthHistogram.cdf :  1 2 3 4 5  6
     *
     * lengths:               10
     * lengthHistogram.data:  10
     * lengthHistogram.cdf :  1 
     */ 

	int SelectRandomValue(T_Data &value) {
		vector<int>::iterator search_it;
        assert(cdf.size() >= 0);
		int randomIndex = RandomInt(cdf[cdf.size()-1]);
		search_it = lower_bound(cdf.begin(), cdf.end(), randomIndex);
		assert(search_it != cdf.end());
		int valueIndex = search_it - cdf.begin();
		value = data[valueIndex];
		return valueIndex;
	}
};
		
#endif
