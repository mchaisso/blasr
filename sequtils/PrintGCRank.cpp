#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"

#include <map>
#include <string>
#include <algorithm>

int MaxHPStretch(Nucleotide *nuc, int length) {
	int p, pn;
	int maxStretch = 0;
	pn = 0;
	while(pn < length) {
		p = pn;
		while(pn < length and toupper(nuc[p]) == toupper(nuc[pn])) { pn ++;}
		if (pn - p > maxStretch) { maxStretch = pn - p;}
	}
	return maxStretch;
}

int MaxSpecHPStretch(Nucleotide *nuc, int length, Nucleotide specNuc) {
	int p, pn;
	pn = 0;
	int maxStretch = 0;
	while(pn < length) {
		while (pn < length and nuc[pn] != specNuc) { pn ++; }
		p = pn;
		while (pn < length and toupper(nuc[p]) == toupper(nuc[pn])) {		pn++;  }
		int stretch = pn - p;
		if (stretch > maxStretch) { maxStretch = stretch; }
	}
	return maxStretch;
}

class RankedIndex {
public:
	float gc;
	int index;
	int maxT, maxA;
	int maxHP;
	int operator<(const RankedIndex &rhs) const {
		return gc > rhs.gc;
	}
	RankedIndex &operator=(const RankedIndex&rhs) {
		gc = rhs.gc; index = rhs.index; maxT = rhs.maxT; maxA = rhs.maxA; maxHP = rhs.maxHP;
		return *this;
	}
};

using namespace std;
int main(int argc, char* argv[]) {
	string infileName;
	int    window;
	if (argc < 3) {
		cout << " usage: printgcrank seq_file window [-mingc gc]" << endl;
		exit(1);
	}

	infileName = argv[1];
	window     = atoi(argv[2]);
	int argi = 3;
	float minGC = 0.0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-mingc") == 0) {
			minGC = atof(argv[++argi]);
		}
		++argi;
	}

	FASTAReader reader;
	reader.Initialize(infileName);
	FASTASequence seq;
	reader.GetNext(seq);

	vector<RankedIndex> rank;
	int p, pi;
	rank.resize(seq.length - window+1);
	for (p = 0; p < seq.length - window + 1; p++) {
		int gc = 0;
		for (pi = p; pi < p + window; pi++) {
			if (toupper(seq.seq[pi]) == 'G' or
					toupper(seq.seq[pi]) == 'C') {
				gc++;
			}
		}
		rank[p].gc = gc*1.0/window;
		rank[p].index = p;
		rank[p].maxT = MaxSpecHPStretch(&seq.seq[pi], window, 'T');
		rank[p].maxA = MaxSpecHPStretch(&seq.seq[pi], window, 'A');
		rank[p].maxHP = MaxHPStretch(&seq.seq[pi], window);
	}
	sort(rank.begin(), rank.end());
	int i; 
	cout << "pct_gc max_t max_a max_hp seq" << endl;
	for (i = 0; i < rank.size(); i++) {
		if (rank[i].gc < minGC) break;
		DNASequence subseq;
		subseq.length = window;
		subseq.seq = &seq.seq[rank[i].index];
		cout << rank[i].gc << " " << rank[i].maxT << " " << rank[i].maxA << " " << rank[i].maxHP << " ";
		subseq.PrintSeq(cout, window+1);
	}
}

