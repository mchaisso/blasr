#include "../common/FASTAReader.h"
#include "../common/DNASequence.h"
#include <string>
#include <map>

using namespace std;
int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "usage: printUngappedLengths in.fasta" << endl;
		exit(1);
	}
	
	string inFileName = argv[1];
	int argi = 2;
	int maxSpan = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-maxSpan") == 0) {
			maxSpan = atoi(argv[++argi]);
		}
		++argi;
	}

	FASTAReader reader;
	reader.Init(inFileName);
	FASTASequence seq;
	int spanIndex = 0;
	cout << "span" <<endl;
	DNALength totSpan = 0;
	while (reader.GetNext(seq)) {
		DNALength p, n;
		DNALength ungappedLength;
		p = 0;
		while (p < seq.length) {
			n = p + 1;
			while (n < seq.length and seq.seq[n] != 'N') {
				n++;
			}
			ungappedLength = n - p;
			if (!maxSpan or ungappedLength < maxSpan) { 
				cout << spanIndex<< " " << ungappedLength << endl;
				totSpan += ungappedLength;
				++spanIndex;
			}
			/*			DNALength bin = ungappedLength/ binWidth;
			if (histogram.find(bin) == histogram.end()) {
				histogram[bin] = 1;
			}
			else {
				histogram[bin]++;
			}
			*/
			while (n < seq.length and seq.seq[n] == 'N') {
				n++;
			}
			p = n;
		}
	}
	cerr << totSpan / spanIndex << endl;
}
