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
	int mapQV;
	AlignmentStats() {
		nMatch = nMismatch = nIns = nDel = 0;
		pctSimilarity = 0.0;
		mapQV = 0;
        score = 0;
	}
	AlignmentStats &Assign(const AlignmentStats &rhs) {
		nMatch = rhs.nMatch;						  		
		nMismatch = rhs.nMismatch;				 
		nIns = rhs.nIns;							  
		nDel = rhs.nDel;							  
		pctSimilarity = rhs.pctSimilarity;
		score = rhs.score;             
		return *this;
	}
 AlignmentStats& operator=(const AlignmentStats &rhs) {
    Assign(rhs);
    return *this;
 }   
  void CopyStats(AlignmentStats rhs) {
    *this = rhs;
  }

};

#endif
