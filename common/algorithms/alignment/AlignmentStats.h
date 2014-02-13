#ifndef ALIGNMENT_STATS_H_
#define ALIGNMENT_STATS_H_


class AlignmentStats {
 public:
	int nMatch;
	int nMismatch;
	int nIns;
	int nDel;
	float pctSimilarity;
	int score;
	AlignmentStats() {
		nMatch = nMismatch = nIns = nDel = 0;
		pctSimilarity = 0.0;
	}
};


#endif
