#ifndef PIECEWISE_MATCH_H_
#define PIECEWISE_MATCH_H_
#include <vector>
using namespace std;
class DirMatch {
 public:
	unsigned int t, q;
	unsigned int s: 1;
	unsigned int v: 1;
	unsigned int l: 30;
	unsigned int len;
	DirMatch() {
		t = q = s = l = 0;
		v = 1;
	}
	bool operator<(const DirMatch &rhs) const {
		if (rhs.t != t) {
			return t < rhs.t;
		}
		else {
			return q < rhs.q;
		}
	}
	DirMatch&operator=(const DirMatch &rhs) {
		t = rhs.t; q = rhs.q; s = rhs.s; v= rhs.v;l=rhs.l; len=rhs.len;
	}

};
		
	

void PiecewiseMatch(  vector<ChainedMatchPos> &forwardMatches,
								 vector<ChainedMatchPos> &reverseMatches,
								 SeqBoundaryFtr<FASTQSequence> &seqBoundary, 
								 WeightedIntervalSet &topIntervals,								 
								 int readLength
								 ) {


	vector<DirMatch> matches;

	matches.resize(forwardMatches.size() + reverseMatches.size());
	if (matches.size() == 0) {
		return;
	}
	int i;
	for (i = 0; i < forwardMatches.size(); i++) {
		matches[i].t = forwardMatches[i].t;
		matches[i].q = forwardMatches[i].q;
		matches[i].s = 0;
		matches[i].l = 0;
		matches[i].len = forwardMatches[i].l;
	}
	int nf=forwardMatches.size();
	for (i = 0; i < reverseMatches.size(); i++) {
		matches[i+nf].t = reverseMatches[i].t;
		matches[i+nf].q = readLength - reverseMatches[i].q;
		matches[i+nf].s = 1;
		matches[i+nf].l = 0;
		matches[i+nf].len = reverseMatches[i].l;
	}
	
	std::sort(matches.begin(), matches.end());


	matches[0].l = 1;

	
	int maxIter = 15;
	int iter;
	vector<int> m(maxIter);
	for (i = 0; i < maxIter; i++) {
		m[i] = (i+1)*3;
	}

	for (iter = 0; iter < maxIter; iter++) {
		for (i = 0; i < matches.size(); i++) {
			matches[i].l = 1;
		}

		i = 0;
		while (i < matches.size()-1) {

			if (matches[i+1].s == matches[i].s and 
					seqBoundary.GetIndex(matches[i].t) == seqBoundary.GetIndex(matches[i+1].t) ) {
				if (matches[i+1].s == 0) {
					if ( matches[i].t < matches[i+1].t and 
							 matches[i].q < matches[i+1].q and 
							 (i+1 == matches.size() - 1 or  
								(matches[i+1].s == matches[i+2].s and 
								 matches[i+1].t < matches[i+2].t and
								 matches[i+1].q < matches[i+2].q))) {
						matches[i+1].l = matches[i].l+1;
					}
					else {
						// restart
						matches[i+1].l = 0;
						matches[i+1].v = 0;
					}
				}
				else {
					if (matches[i].t < matches[i+1].t and 
							matches[i].q > matches[i+1].q and
							(i+1 == matches.size()-1 or 
							 (matches[i+1].s == matches[i+2].s and 
								matches[i+1].t < matches[i+2].t and
								matches[i+1].q > matches[i+2].q))) {
						matches[i+1].l = matches[i].l+1;
					}
					else {
						matches[i+1].l = 0;
						matches[i+1].v = 0;
					}
				}
			}
			else {
				if (i < matches.size()-1) {
					matches[i+1].l = 0;
					matches[i+1].v = 0;
				}
			}
			i = i+1;
		}

		i = 0;
		int start;
		int n;
		int nr = 0;
		int j;
		while (i < matches.size()) {
			while (i < matches.size() and matches[i].v == 0) { i++; }
			start = i;
			j = i + 1;
			n = 1;
			while (j < matches.size() and matches[j].l != 0) {
				while (j < matches.size() and matches[j].v == 0) { j++; }
				j++;
				n++;
			}
			if (n < m[iter]) {
				int t;
				for (t = start; t < j; t++) {
					matches[t].v = 0;
				}
				nr += n;
			}
			i = j;
		}
		i = j = 0;
		for (i = 0; i < matches.size(); i++) {
			if (matches[i].v != 0) { assert(j < matches.size());matches[j] = matches[i]; j++;}

		}
		if (j > 0) {
			matches.resize(j);
		}
		else { 
			matches.clear();
			break;
		}
	}

	int start, end;
	i = 0;
	if (matches.size() == 0) {
		return;
	}
	while (i < matches.size()-1) {
		start = i;
		while (i < matches.size() -1 and matches[i+1].l > matches[i].l and 
					 abs((int)(matches[i+1].t - matches[i].t)) < 30000 and 
					 abs((int)(matches[i+1].q - matches[i].q)) < 30000
					 ) { i++; }
		end = i;
		i++;
		int qStart  = matches[start].q;
		int qEnd    = matches[end].q + matches[end].len;

		if (matches[start].s == 1) {
			qEnd = readLength - matches[end].q;
			qStart = readLength - matches[start].q;
		}


		vector<ChainedMatchPos> emptyMatches;

		WeightedInterval weightedInterval(end-start, (unsigned int) end - start, (unsigned int) matches[end].q - matches[start].q,
																			matches[start].t,  matches[end].t + matches[end].len,
																			(int) matches[start].s, 1.0/(end-start),
																			qStart, qEnd, readLength, emptyMatches);
		topIntervals.insert(weightedInterval);		
	}
	
}
	


#endif
