#ifndef VARIABLE_LENGTH_SDP_FRAGMENT_H_
#define VARIABLE_LENGTH_SDP_FRAGMENT_H_
#include "SDPFragment.h"

class ChainedFragment : public Fragment {
	int score;
	ChainedFragment *chainPrev;
 public:
    ChainedFragment():Fragment() {
        score = 0;
        chainPrev = NULL; 
    }
	
	int SetScore(int _score) {
		return (score = _score);
	}

	int GetScore() {
		return score;
	}
	
	ChainedFragment *SetChainPrev(ChainedFragment *_chainPrev) {
		return (chainPrev = _chainPrev);
	}

	ChainedFragment *GetChainPrev() {
		return chainPrev;
	}
	
	int LessThan(const ChainedFragment &f) const {
		//
		// Order fragments by endpoint.
		//
		if (x == f.GetX()) {
			return (y < f.GetY());
		}
		else {
			return x < f.GetX();
		}
	}

	int operator<(const ChainedFragment &f) const {
		// 
		// Sort fragments by x then y.
		//
		return LessThan(f);
	}
};


#endif
