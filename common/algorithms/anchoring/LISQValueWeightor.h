#ifndef ALGORITHMS_ANCHORING_LIS_QVALUE_WEIGHTOR_H_
#define ALGORITHMS_ANCHORING_LIS_QVALUE_WEIGHTOR_H_

#include "../../qvs/QualityValue.h"
#include "../../DNASequence.h"
#include "../../FASTQSequence.h"

template<typename T_MatchList, typename T_Sequence>
class LISQValueWeightor {
 public:
	T_Sequence *seq;
	float operator()(const T_MatchList &matchList) {
		float totalQ;
		DNALength  nBases;
		VectorIndex i;
		totalQ = 0.0;
		nBases = 0;
		for (i = 0; i < matchList.size(); i++) {
			DNALength mp;
			for (mp = matchList[i].q; mp < matchList[i].q + matchList[i].w; mp++) {
				totalQ += (*seq).qual[mp];
			}
			nBases += matchList[i].w;
		}
		if (nBases > 0) {
			return totalQ / nBases;
		}
	}
};
			


#endif
