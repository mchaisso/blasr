#include "algorithms/sorting/RadixSort.h"
#include "FASTAReader.h"
#include "FASTASequence.h"

#include <string>

int main(int argc, char* argv[]) {

	if (argc < 3) {
		cout << "usage: testRadixSort inFile numDigits [-suffix] [-words] [-stable]" << endl;
	}

	string inFile = argv[1];
	unsigned int numDigits = atoi(argv[2]);
	int argi = 3;
	int sortType = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-suffix") == 0) {
			sortType = 0;
		}
		else if (strcmp(argv[argi], "-words") == 0) {
			sortType = 1;
		}
	else if (strcmp(argv[argi], "-stable") == 0) {
			sortType = 2;
		}
		++argi;
	}
	FASTAReader reader;
	reader.Init(inFile);
	FASTASequence seq;
	reader.GetNext(seq);
	seq.ToTwoBit();
	if (seq.length < numDigits) {
		cout << "no keys" << endl;
		exit(0);
	}
	int *v = new int[seq.length];
	DNALength i;

	for (i = 0; i < seq.length; i++) {	v[i] = i;}
	Nucleotide *key = new Nucleotide[numDigits + 1];
	/*
	cout << "before in-place: " << endl;

	for (i = 0; i < seq.length; i++) {
		int j;
		int wordLength = numDigits;
		if (v[i] + numDigits > seq.length) {
			wordLength = seq.length - v[i];
		}
		key[wordLength] = '\0';
		for (j = 0; j < wordLength; j++ ) {key[j] = TwoBitToAscii[seq.seq[v[i] + j]];}
		cout << key << endl;
	}
	*/	
	key = new Nucleotide[numDigits + 1];
	if (sortType == 0) {
		RadixSuffixSortLR(v, seq.seq, seq.length, numDigits, 4);

		cout << "result of in-place: " << endl;

		for (i = 0; i < seq.length; i++) {
			unsigned int j;
			unsigned int wordLength = numDigits;
			if (v[i] + numDigits > seq.length) {
				wordLength = seq.length - v[i];
			}
			key[wordLength] = '\0';
			for (j = 0; j < wordLength; j++ ) {key[j] = TwoBitToAscii[seq.seq[v[i] + j]];}
			cout << key << endl;
		}
	}
	else if (sortType == 1 or sortType == 2) {

		unsigned int numKeys = seq.length - numDigits + 1;
		if (sortType == 1) {
			RadixSortLR(v, seq.seq, numKeys, numDigits, 4);
		}
		else {
			int *v2 = new int[seq.length];
			RadixSortRL(v, v2, seq.seq, numKeys, numDigits, 4);
		}
		for (i = 0; i < numKeys; i++) {
			int j;
			int wordLength = numDigits;
			if (v[i] + numDigits > seq.length) {
				wordLength = seq.length - v[i];
			}
			key[wordLength] = '\0';
			for (j = 0; j < wordLength; j++ ) {key[j] = TwoBitToAscii[seq.seq[v[i] + j]];}
			cout << key << endl;
		}
	}


}
	
