#include <string>
#include <string.h>
#include <algorithm>
#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"

class KCompareFunctor {
	
public:
	Nucleotide *seq;
	int k;
	int operator()(unsigned int a, unsigned int b) {
		return strncmp((const char*) &seq[a], (const char*) &seq[b], k) < 0;
	}
};


int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: countKmer in.fa k" << endl;
		exit(1);
	}
	string sequenceFileName = argv[1];
	int k = atoi(argv[2]);

	unsigned int *indexArray;
	
	FASTAReader reader;
	reader.Init(sequenceFileName);
	FASTASequence seq;
	reader.GetNext(seq);
	seq.ToUpper();
	unsigned int nKmers = seq.length - k + 1;
	indexArray = new unsigned int[nKmers];
	unsigned int i;
	for (i = 0; i < nKmers; i++ ) {
		indexArray[i] = i;
	}
	KCompareFunctor comp;
	comp.k = k;
	comp.seq = seq.seq;
	std::sort(indexArray, indexArray + nKmers, comp);
	unsigned int nUnique = 0;
	for (i = 0; i < nKmers-1; ) {
	  int j = i + 1;
	  while (j < nKmers and strncmp((const char*) &seq.seq[indexArray[i]], (const char*) &seq.seq[indexArray[j]], k) == 0) j++;
	  if (j == i + 1) {
			nUnique++;
		}
	  i = j;
	}
	cout << nUnique << " " << nKmers << endl;
	return 0;
}
