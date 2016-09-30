#ifndef MATCH_POS_H_
#define MATCH_POS_H_

#include <vector>
#include <algorithm>
#include "DNASequence.h"

using namespace std;


class MatchPos {
 public:
	DNALength t, q;
	MatchWeight w;
	DNALength l;
	MatchPos(DNALength pt, DNALength pq, DNALength pl) {
		t = pt;
		q = pq;
		l = pl;
		w = 0;
	}
	MatchPos (const MatchPos &rhs) {
		(*this) = rhs;
	}
	MatchWeight Size() {
		return l;
	}
	
	MatchPos() {
		t = q = -1;
		l = 0;
		w = 0;
	}

	MatchPos& operator=(const MatchPos &rhs) {
		t = rhs.t; q = rhs.q; w = rhs.w;
		l = rhs.l;
		return *this;
	}
	
	DNALength GetLength() const {
		return l;
	}
	DNALength GetXLength() const {
		return l;
	}
	DNALength GetYLength() const {
		return l;
	}

	DNALength GetX() const {
		return q;
	}
	DNALength GetY() const {
		return t;
	}
	UInt GetT() {
		return t;
	}

	UInt GetQ() {
		return (UInt) q;
	}

	UInt GetTW() {
		return w;
	}
	
	UInt GetQW() {
		return w;
	}
	
	UInt GetW() {
		return w;
	}
  friend ostream& operator<<(ostream & out, MatchPos &p) {
    out << p.q << "\t" << p.t <<"\t"<< p.l;
    return out;
  }
};


class ChainedMatchPos : public MatchPos {
 protected:
	int score;
	ChainedMatchPos *chainPrev;
 public:
 ChainedMatchPos(DNALength pt, DNALength pq, DNALength pl) : MatchPos(pt, pq, pl) {score = 0; chainPrev = NULL;}
  ChainedMatchPos() : MatchPos() {
		score = 0;
		chainPrev = NULL;
	}
	ChainedMatchPos(const ChainedMatchPos &rhs) {
		(*this) = rhs;
	}
	int GetScore() {
		return score;
	}
		
	int SetScore(int _score) {
		return (score = _score);
	}
	ChainedMatchPos* SetChainPrev(ChainedMatchPos *_chainPrev) {
		return (chainPrev = _chainPrev);
	}
	ChainedMatchPos* GetChainPrev() {
		return chainPrev;
	}
	ChainedMatchPos &operator=(const ChainedMatchPos &rhs) {
		((MatchPos&)(*this)) = ((const MatchPos&)rhs);
		return *this;
	}
  friend ostream& operator<<(ostream & out, ChainedMatchPos &p) {
    out << p.q << "\t" << p.t <<"\t"<< p.l;
    return out;
  }
};

template<typename T_MatchPos>
class CompareMatchPos {
 public:
	int operator()(const T_MatchPos &lhs, const T_MatchPos &rhs) const {
		if (lhs.t < rhs.t) 
			return 1;
		else if (lhs.t > rhs.t)
			return 0;
		else { 
			return lhs.q < rhs.q;
		}
	}
};


typedef vector<MatchPos> MatchPosList;

template<typename T_MatchPos>
class CompareMatchPosByWeight {
 public:
	int operator()(const T_MatchPos &a, const T_MatchPos &b) const {
		return a.l < b.l;
	}
};


template<typename T_MatchPos>
class CompareMatchPosIndexByWeight {
 public:
	vector<T_MatchPos> *list;
	int operator()(const int i, const int j) const {
		return ((*list)[i].l > (*list)[j].l);
	}
};


template<typename T_MatchPos>
class CompareMatchPosIndexByTextPos {
 public:
	vector<T_MatchPos> *list;
	int operator()(const int i, const int j) const {
		return (*list)[i].t < (*list)[j].t;
	}
};


template<typename T_MatchPos>
void SortMatchPosList(vector<T_MatchPos> &mpl) {
	std::sort(mpl.begin(), mpl.end(), CompareMatchPos<T_MatchPos>());
}

template<typename T_MatchPos>
void SortMatchPosListByWeight(vector<T_MatchPos> &mpl) {
	std::sort(mpl.begin(), mpl.end(), CompareMatchPosByWeight<T_MatchPos>());
}

template<typename T_MatchPos>
void SortMatchPosIndexListByWeight(vector<T_MatchPos> &mpl, vector<int> &indices) {
	CompareMatchPosIndexByWeight<T_MatchPos> cmp;
	cmp.list = &mpl;
	std::sort(indices.begin(), indices.end(), cmp);
}

template<typename T_MatchPos>
void SortMatchPosIndexListByTextPos(vector<T_MatchPos> &mpl, vector<int> &indices) {
	CompareMatchPosIndexByTextPos<T_MatchPos> cmp;
	cmp.list = &mpl;
	std::sort(indices.begin(), indices.end(), cmp);
};

#endif
