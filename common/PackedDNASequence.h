#ifndef PACKED_DNA_SEQUENCE_H_
#define PACKED_DNA_SEQUENCE_H_
#include <limits.h>
#include <vector>

#include "defs.h"
#include "DNASequence.h"
#include "NucConversion.h"

typedef unsigned int t_word;
class PackedDNASequence {
 public:
	vector<t_word> pSeq;
	int nucPerWord;
	int Pack(DNASequence &seq) {
		int pos = 0;
		int numPacked;
		nucPerWord = __WORDSIZE / dna.bitsPerNuc;
		pSeq.resize(seq.length / nucPerWord);
		
		while (pos < seq.length) {
			numPacked = PackWord

	}

	int PackWord(DNASequence &dna, int length, t_word word) {
		int i;

		word = 0;
		int nucInWord = MIN(nucPerWord, length);
		if (nucInWord == 0) {
			return 0;
		}
		for (i = 0; i < nucInWord - 1; i++) {
			word += dna.seq[i];
			word = word << ((t_word) dna.bitsPerNuc);
		}
		word += dna.seq[i];
		return nucInWord;
	}
};




#endif
