#include "datastructures/suffixarray/SuffixArray.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include <string>
#include <vector>

using namespace std;
int main(int argc, char* argv[]) {
	string genomeFileName;
	string suffixArrayFileName;
	if (argc < 4) {
		cout << "Usage: countKMers genome suffixArray k [k2 k3 k4...]" << endl;
		exit(1);
	}
	genomeFileName = argv[1];
	suffixArrayFileName = argv[2];
	int argi = 3;
	vector<DNALength> k;
	while (argi < argc) {
		k.push_back(atoi(argv[argi]));
		argi++;
	}

	// Get the ref sequence.
	FASTAReader reader;
	reader.Init(genomeFileName);
	FASTASequence seq;
	reader.GetNext(seq);
	seq.ToUpper();
	// Get the suffix array.
	DNASuffixArray sarray;
	sarray.Read(suffixArrayFileName);
	
	int ki;
	for (ki = 0; ki < k.size(); ki++) {

		DNALength i;
		DNALength numUnique = 0;
		for (i = 0; i < seq.length - k[ki] - 1; ) {
			DNALength j = i + 1;
			while (j < seq.length - k[ki] and 
						 seq.length - sarray.index[i] >= k[ki] and
						 seq.length - sarray.index[j] >= k[ki] and 
						 strncmp((const char*) &seq.seq[sarray.index[i]], (const char*) &seq.seq[sarray.index[j]], k[ki]) == 0) {
				j++;
			}
			if (j == i + 1) { 
				++numUnique;
			}
			i = j;
		}
		cout << k[ki] << " " << numUnique << " " << sarray.length << endl;
	}
}
