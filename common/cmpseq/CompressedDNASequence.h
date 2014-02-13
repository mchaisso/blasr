#ifndef COMPRESSED_DNA_SEQUENCE_H_
#define COMPRESSED_DNA_SEQUENCE_H_

#include "DNASequence.h"
#include "FASTASequence.h"
#include "algorithms/compare/Compare4BitCompressed.h"
#include "defs.h"
#include <string.h>
#include "qvs/QualityValue.h"

typedef unsigned char CompressedNucleotide;

class CompressedDNASequence: public DNASequence {
	static const unsigned char MaskCount = 0xf;
	static const unsigned char MaskNuc   = 0xf0;
	static const unsigned char ShiftCount = 4;

 public:
	char *title;
	int titleLength;
	//
	// This is just a placeholder for now.  
	// No extra data here, just the ability to decompress.  Right now 
	// the utilities for the compressed dna sequences
	// are in CompressedSeqUtils.h, which could move here later.
	//
	QualityValue *qual;

	CompressedDNASequence() {
		const char t[] = "Compressed sequence\0";
		titleLength = strlen(t);
		title = new char[titleLength+1];
		strcpy(title, t);
		title[titleLength] = '\0';
	}
		
	void MakeRC(CompressedDNASequence &rc) {
		rc.Allocate(length);
		DNALength i;
		for (i = 0; i < length; i++) {
			rc.seq[length - i - 1] = ReverseComplementNuc[ThreeBit[seq[i] & MaskCount]];
			rc.seq[length - i - 1] += (seq[i] & MaskNuc); 
		}
		memcpy(rc.title, title, titleLength);
		rc.titleLength = titleLength;
	}
	
	Nucleotide operator[](DNALength i) {
		return GetNuc(i);
	}

	Nucleotide GetNuc(DNALength i) {
		return (seq[i] & MaskCount);
	}

	unsigned char GetCount(DNALength i) {
		return seq[i] >> ShiftCount;
	}

	char *GetName() {
		return (char*) title;
	}
	
	void Copy(FASTASequence &rhs) {
		seq = new CompressedNucleotide[rhs.length];
		memcpy(seq, rhs.seq, rhs.length);
		length = rhs.length;
		if (title != NULL) {
			delete[] title;
		}
		title = new char[rhs.titleLength+1];
		memcpy(title, rhs.title, rhs.titleLength);
		titleLength = rhs.titleLength;
		title[titleLength] = '\0';
	}
	float GetAverageQuality() {
		return 0.0;
	}

	void SortHomopolymerQualities() {
		cout << "qualities are not implemented for compressed sequences." << endl;
		assert(0);
	}
};

	


#endif
