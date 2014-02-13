#include "datastructures/suffixarray/SuffixArray.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include <string>


using namespace std;
int main(int argc, char* argv[]) {
	string genomeFileName;
	string suffixArrayFileName;
	if (argc < 4) {
		cout << "Usage: printWordCount genome suffixArray k [k2 k3 k4...]" << endl;
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
  //	reader.GetNext(seq);
  reader.ReadAllSequencesIntoOne(seq);
	seq.ToUpper();
	// Get the suffix array.
	DNASuffixArray sarray;
	sarray.Read(suffixArrayFileName);
	
	int ki;
  char *word;
  cout << "wordlen word nword" << endl;
	for (ki = 0; ki < k.size(); ki++) {
    word = new char[k[ki]+1];
    word[k[ki]] = '\0';
		DNALength i;
		DNALength numUnique = 0;
		for (i = 0; i < seq.length - k[ki] - 1; ) {
			DNALength j = i + 1;
      bool seqAtN = false;
      int si;
      for(si = 0; si < k[ki]; si++) {
        if (seq.seq[sarray.index[i] + si] == 'N') {
          seqAtN = true;
          break;
        }
      }
      if (seqAtN) {
        i++;
        continue;
      }
			while (j < seq.length - k[ki] and 
						 seq.length - sarray.index[i] >= k[ki] and
						 seq.length - sarray.index[j] >= k[ki] and 
						 strncmp((const char*) &seq.seq[sarray.index[i]], (const char*) &seq.seq[sarray.index[j]], k[ki]) == 0) {
				j++;
			}
      if (seq.length - sarray.index[i] >= k[ki]) {
        for(si = 0; si < k[ki]; si++) {
          word[si] = seq.seq[sarray.index[i]+si];
        }
        cout << k[ki] << " " << word << " " << j - i + 1 << endl;
        if (j == i + 1) { 
          ++numUnique;
        }
      }
			i = j;
		}
	}
}
