#ifndef DATASTRUCTURES_BWT_BWT_H_
#define DATASTRUCTURES_BWT_BWT_H_

#include <iostream>
#include <fstream>
#include "Occ.h"
#include "Pos.h"
#include "../../datastructures/suffixarray/SuffixArray.h"
#include "../../datastructures/sequence/PackedDNASequence.h"
#include "../../FASTASequence.h"

using namespace std;

/*
 * Define an Occurrence table appropriate for Gb sized genomes.
 * Probably everything will end up using this.
 */
typedef Occ <PackedDNASequence, unsigned int, unsigned short> GbOcc;

/*
 * Define an Occurrence table appropriate for Mb sized genomes.
 */
typedef Occ <PackedDNASequence, unsigned int, unsigned char> MbOcc;


class SingleStoragePolicy {
 public:
	DNALength *spp, *epp;
	void Store(DNALength sp, DNALength ep) {
		*spp = sp;
		*epp = ep;
	}
};

class VectorStoragePolicy {
 public:
	vector<DNALength> *spvp, *epvp;
	void Store(DNALength sp, DNALength ep) {
		spvp->push_back(sp);
		epvp->push_back(ep);
	}
};

template<typename T_BWT_Sequence, typename T_DNASequence>
class Bwt {
 public:
	T_BWT_Sequence bwtSequence;
	GbOcc occ;
	Pos<T_BWT_Sequence>   pos;
	static const int CharCountSize = 7;
	int useDebugData;
	vector<DNALength> saCopy;
	DNALength charCount[CharCountSize];
	DNALength firstCharPos;
	

	void Write(string outName) {
		ofstream bwtOut;
		CrucialOpen(outName, bwtOut, std::ios::binary|std::ios::out);
		Write(bwtOut);
	}

	void PrintBWTString(ostream &out) {
		DNALength p;
		for (p = 0; p < bwtSequence.length; p++) {
			out << (char) ThreeBitToAscii[bwtSequence[p]];
			if (p % 50 == 49) out << endl;
		}
		if(p % 50 != 0) out << endl;
	}

	void Write(ostream &bwtOut) {
		bwtSequence.Write(bwtOut);
		bwtOut.write((char*)charCount, sizeof(DNALength)*CharCountSize);
		bwtOut.write((char*)&firstCharPos, sizeof(DNALength));
		bwtOut.write((char*)&useDebugData, sizeof(useDebugData));
		if (useDebugData) {
			bwtOut.write((char*)&saCopy[0], (bwtSequence.length-1) * sizeof(DNALength));
		}
		occ.Write(bwtOut);
		pos.Write(bwtOut);
	}

	int Read(string inName) {
		ifstream bwtIn;
		DNALength seqStorageSize;
		CrucialOpen(inName, bwtIn, std::ios::binary|std::ios::in);
		bwtSequence.Read(bwtIn);
		bwtIn.read((char*)charCount, sizeof(DNALength)*CharCountSize);
		bwtIn.read((char*)&firstCharPos, sizeof(DNALength));
		bwtIn.read((char*)&useDebugData, sizeof(useDebugData));
		if (useDebugData) {
			saCopy.resize(bwtSequence.length-1);
			bwtIn.read((char*)&saCopy[0], (bwtSequence.length-1) * sizeof(DNALength));
		}
		occ.Read(bwtIn, useDebugData);
		pos.Read(bwtIn);
		occ.InitializeBWT(bwtSequence);
		return 1;
	}

	void Print(ofstream &out) {
		bwtSequence.Print(out);
	}

	DNALength LFBacktrack(DNALength bwtPos) {
		Nucleotide curNuc = bwtSequence.Get(bwtPos);
		assert(curNuc < 5);
		DNALength bwtPrevPos = charCount[curNuc] + occ.Count(curNuc, bwtPos) - 1;
		return bwtPrevPos;
	}

	DNALength Locate(DNALength bwtPos) {
		DNALength seqPos;
		DNALength offset = 0;
		while (1) {
			if (pos.Lookup(bwtPos, seqPos)) {
				break;
			}
			else {
				DNALength bwtPrevPos;
				bwtPrevPos = LFBacktrack(bwtPos);
				if (useDebugData) {
					assert(saCopy[bwtPos-1] - 1 == saCopy[bwtPrevPos-1]);
				}
				bwtPos = bwtPrevPos;
				assert(bwtPos <= bwtSequence.length);
				/*
				 * Boundary condition at the beginning of the bwt string.
				 */
				if (bwtPos == firstCharPos) {
					seqPos = 1;
					break;
				}
			}
			++offset;
		}
		return seqPos + offset;
	}
	
	DNALength Locate(DNALength sp, DNALength ep, vector<DNALength> &positions, int maxCount = 0) {
		DNALength bwtPos;
		DNALength seqPos;
		if (sp <= ep and (maxCount == 0 or ep - sp < maxCount)) {
			for (bwtPos = sp; bwtPos <= ep; bwtPos++) {
				if ((seqPos = Locate(bwtPos))) {
					positions.push_back(seqPos);
				}
			}
		}
		return ep - sp + 1;
	}

	DNALength Locate(T_DNASequence &seq, vector<DNALength> &positions, int maxCount =0) {
		DNALength ep, sp;
		Count(seq, sp, ep);
		return Locate(sp, ep, positions);
	}

	DNALength GetNumCharsLessThan(Nucleotide nuc) {
		return charCount[nuc];
	}


	void InitializeDNACharacterCount() {
		//
		// All counts start at 1 due to implicit encoding of $ character,
		// where $ is less than all other chars.
		//
		fill(charCount, &charCount[6], 0);

		DNALength p;
		Nucleotide nuc;
		for (p = 0; p < bwtSequence.length; p++) {
			nuc = bwtSequence[p];
				/*
				 * Intentionally omit break commands so that charCount[4]
				 * contains the counts of all characters 4 and below and so on.
				 */
			switch(nuc) {
		  case 5:
				//
				// 5 is out of order here because the '$' is stored after ACGTN
				// since it is a nonstandard character.
				//
				charCount[0]++;
			case 0: //A
				charCount[1]++;
			case 1: //C
				charCount[2]++;
			case 2: //G
				charCount[3]++;
			case 3: //T
				charCount[4]++;
			case 4: //N
				charCount[5]++;
			}
		}
		// sum
		charCount[6] = bwtSequence.length;
	}

	template<typename T_PStoragePolicy>
	int Count(T_DNASequence &seq, T_PStoragePolicy &StoragePolicy) {
		
		/*
		 * Implement algorithm count directly from the FM-Index paper(s --
		 * it's shown many times).
		 */
		DNALength p = seq.length-1;
		DNALength sp, ep;
		Nucleotide c;
		int i;

		//
		// Original forumlation is using count offsets starting at 1.
		//
		Nucleotide tbn = ThreeBit[seq[p]];
		sp = charCount[tbn]; // = +1 (from paper) - 1 (0
		                     // offset not in paper).
		ep = charCount[tbn +1] - 1;
		StoragePolicy.Store(sp,ep);
		while (sp <= ep and p > 0) {
			c  = ThreeBit[seq[p-1]];    
			int  cc = charCount[c];
			sp = cc + occ.Count(c,sp-1) + 1 - 1;
			ep = cc + occ.Count(c,ep) - 1;
			StoragePolicy.Store(sp,ep);
			p--;
		}
		return ep - sp + 1;
	}

	int Count(T_DNASequence &seq, DNALength &sp, DNALength &ep) {
		/*
		 * Implement algorithm count directly from the FM-Index paper(s --
		 * it's shown many times).
		 */
		SingleStoragePolicy storagePolicy;
		storagePolicy.spp = &sp;
		storagePolicy.epp = &ep;
		return Count(seq, storagePolicy);
	}
	
	int Count(T_DNASequence &seq, vector<DNALength> &spv, vector<DNALength> &epv) {
		VectorStoragePolicy storagePolicy;
		storagePolicy.spvp = &spv;
		storagePolicy.epvp = &epv;
		return Count(seq, storagePolicy);
	}

	int Count(T_DNASequence &seq) {
		DNALength ep, sp;
		return Count(seq, sp, ep);
	}

	void InitializeBWTStringFromSuffixArray(T_DNASequence &origSeq, DNALength saIndex[]) {
		// extra +1 is for $.
		bwtSequence.Allocate(origSeq.length+1);
		if (useDebugData) {
			saCopy.resize(origSeq.length);
		}

		DNALength p;

		if (origSeq.length == 0) {
			//
			// No work to do, but even the null string has the sentinal
			// appended to it.
			//
			bwtSequence.Set(0,ThreeBit[(int)'$']);
			return;
		}
		
		//
		// By convention, bwt[0] = T[len(T)-1] because T[len(T)] == '$',
		// the lexicographic least character in the alphabet.
		//
		bwtSequence.Set(0, ThreeBit[origSeq[origSeq.length-1]]);
		
		for (p = 1; p < origSeq.length+1; p++) {
			if (useDebugData) {
				saCopy[p-1] = saIndex[p-1];
			}
			if (saIndex[p-1] > 0) {
				assert(ThreeBit[origSeq[saIndex[p-1]-1]] != 255);
				bwtSequence.Set(p, ThreeBit[origSeq[saIndex[p-1]-1]]);
			}
			else {
				//
				// The 0'th suffix corresponds to the one ending in the
				// sentinal '$'.  Since this is explicitly encoded, we can
				// store a value in the bwt for '$'.
				//
				firstCharPos = p;
				bwtSequence.Set(p,ThreeBit[(int)'$']);
			}
		}
	}										

	void InitializeFromSuffixArray(T_DNASequence &dnaSeq, DNALength saIndex[], int buildDebug=0) {
		useDebugData = buildDebug;
		InitializeBWTStringFromSuffixArray(dnaSeq, saIndex);
		InitializeDNACharacterCount();

		// sequence, major, minor bin sizes.
		occ.Initialize(bwtSequence, 4096, 64, buildDebug);
		pos.InitializeFromSuffixArray(saIndex, dnaSeq.length);
	}
};

typedef	Bwt<PackedDNASequence, FASTASequence> BWT;


#endif
