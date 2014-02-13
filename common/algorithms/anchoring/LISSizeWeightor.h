#ifndef LIS_SIZE_WEIGHTOR_H_
#define LIS_SIZE_WEIGHTOR_H_

#include "../../Types.h"

template<typename T_MatchList>
class LISSizeWeightor {
 public:
	MatchWeight operator()(T_MatchList &matchList) {
		MatchWeight size = 0;
		VectorIndex i;
		for (i = 0; i < matchList.size(); i++) {
			//			size += matchList[i].GetWeight();
			size += matchList[i].GetLength();
		}
		return size;
	}
};


#endif
