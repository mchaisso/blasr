#ifndef DATASTRUCTURES_BWT_OCC_H_
#define DATASTRUCTURES_BWT_OCC_H_

#include "../../DNASequence.h"
#include "../../NucConversion.h"
#include "../../utils.h"
#include "../matrix/Matrix.h"
#include "../matrix/FlatMatrix.h"

template<typename T_BWTSequence, typename T_Major, typename T_Minor>
class Occ {
 public:
	int majorBinSize;
	int minorBinSize;
	int hasDebugInformation;
	FlatMatrix2D<T_Major> major;
	FlatMatrix2D<T_Minor> minor;
	FlatMatrix2D<T_Major> full;
	static const unsigned int AlphabetSize = 5;
	T_BWTSequence *bwtSeqRef;
	DNALength numMajorBins, numMinorBins;
	void PrintBins(ostream &out) {
		out << "numMajor: " << numMajorBins << " numMinor: " << numMinorBins << endl;
		DNALength ma, mi, mii;
		mi = 0;
		int i;
		for (ma = 0; ma < numMajorBins; ma++) {
			out << "MAJOR: ";
			for (i = 0; i < 5; i++ ){ out << major[ma][i] << " "; } out << endl;
			for (mii = 0; mii < majorBinSize / minorBinSize  && mi < numMinorBins; mii++, mi++) {
				out << "       ";
				for (i = 0;i <5; i++ ){
					out << minor[mi][i] << " ";
				}
				out << endl;
			}
		}
	}
	
	void  InitializeBWT(T_BWTSequence &bwtSeq) {
		bwtSeqRef    = &bwtSeq;
	}

	int Initialize(T_BWTSequence &bwtSeq,
								 int _majorBinSize=4096,
								 int _minorBinSize=64,
								 int _hasDebugInformation=0) {
		//
		// This reference is used when counting nucleotides. It assumes
		// the sequence does not change between initialization and
		// subsequent calls to count.
		//
		bwtSeqRef    = &bwtSeq;

		majorBinSize = _majorBinSize;
		minorBinSize = _minorBinSize;
		hasDebugInformation = _hasDebugInformation;
		InitializeMajorBins(bwtSeq);
		InitializeMinorBins(bwtSeq);
		if (hasDebugInformation) {
			InitializeTestBins(bwtSeq);
		}
	}

	void InitializeMajorBins(T_BWTSequence &bwtSeq) {
		numMajorBins = CeilOfFraction(bwtSeq.length, (DNALength) majorBinSize);
		major.Allocate(numMajorBins, AlphabetSize);
		vector<DNALength> runningTotal;
		runningTotal.resize(AlphabetSize);
		fill(runningTotal.begin(), runningTotal.end(), 0);
		fill(&major.matrix[0], &major.matrix[numMajorBins*AlphabetSize], 0);
		DNALength p;
		DNALength binIndex = 0;
		for (p = 0; p < bwtSeq.length; p++) {
			Nucleotide nuc = ThreeBit[bwtSeq[p]];
			// only handle ACTGN, $==6, so skip counting that.
			if (nuc > AlphabetSize) continue;
			if (p % majorBinSize == 0) { //majorBinSize-1) {
				//				cout << "storing at " << p<< " " << binIndex << endl;
				int n;
				for (n = 0; n < AlphabetSize; n++ ) {
					major[binIndex][n] = runningTotal[n];
				}
				binIndex++;
			}
			runningTotal[nuc]++;
		}
	}

	void InitializeTestBins(T_BWTSequence &bwtSeq) {
		full.Allocate(bwtSeq.length, AlphabetSize);
		fill(full.matrix, &full.matrix[bwtSeq.length * AlphabetSize],0);
		DNALength p;
		int n;
		for (p = 0; p < bwtSeq.length; p++) {
			Nucleotide nuc = ThreeBit[bwtSeq[p]];
			if (nuc > AlphabetSize) {
				for (n = 0; n < AlphabetSize; n++ ) {
					full[p][n] = full[p-1][n];
				}
			}
			else {
				full[p][nuc]++;
				if (p > 0) {
					for (n = 0; n < AlphabetSize; n++ ) {
						full[p][n] = full[p-1][n] + full[p][n];
					}
				}
			}
		}
	}

	void InitializeMinorBins(T_BWTSequence &bwtSeq) {
		numMinorBins = CeilOfFraction(bwtSeq.length, (DNALength) minorBinSize);
		minor.Allocate(numMinorBins, AlphabetSize);

		vector<DNALength> majorRunningTotal;
		majorRunningTotal.resize(AlphabetSize);		
		fill(majorRunningTotal.begin(), majorRunningTotal.end(), 0);
		fill(&minor.matrix[0], &minor.matrix[numMinorBins*AlphabetSize], 0);
		
		DNALength p;
		DNALength minorBinIndex = 0;
		for (p = 0; p < bwtSeq.length; p++ ){
			Nucleotide nuc = ThreeBit[bwtSeq[p]];
			if (nuc > AlphabetSize) continue;
			//
			//  The minor bins are running totals inside each major
			//  bin. When the count hits a bin offset, reset the bin
			//  counter. 
			//  
			if (p % majorBinSize == 0) {
				fill(majorRunningTotal.begin(), majorRunningTotal.end(), 0);				
			}
			if (p % minorBinSize == 0) {
				int n;
				for (n = 0; n < AlphabetSize; n++ ) {
					minor[minorBinIndex][n] = majorRunningTotal[n];
				}
				minorBinIndex++;
			}
			majorRunningTotal[nuc]++;
		}
	}

	int Count(Nucleotide nuc, DNALength p ) {
		DNALength majorIndex = p / majorBinSize;
		DNALength minorIndex = p / minorBinSize;
		DNALength lastBinnedIndex = minorBinSize * (p / minorBinSize);
		//
		// This should be sort of O(1), since the last expression should
		// be made of bit operations that are fast.
		//
		Nucleotide smallNuc = ThreeBit[nuc];
		//assert(smallNuc < 5);
		DNALength nocc = major[majorIndex][smallNuc] + minor[minorIndex][smallNuc] + bwtSeqRef->CountNuc(lastBinnedIndex, p+1, nuc);
		//		assert(full.matrix == NULL or full[p][smallNuc] == nocc);
		return nocc;
	}

	void Write(ostream &out) {
		out.write((char*) &majorBinSize, sizeof(majorBinSize));
		out.write((char*) &minorBinSize, sizeof(minorBinSize));
		
		out.write((char*) &numMajorBins, sizeof(numMajorBins));
		if (numMajorBins > 0) {
			out.write((char*) major[0], sizeof(T_Major) * numMajorBins*AlphabetSize);
		}
		
		out.write((char*) &numMinorBins, sizeof(numMinorBins));
		if (numMinorBins > 0) {
			out.write((char*) minor[0], sizeof(T_Minor) * numMinorBins*AlphabetSize);
		}
		if (hasDebugInformation) {
			DNALength bwtSeqLength = bwtSeqRef->length;
			out.write((char*)&bwtSeqLength, sizeof(bwtSeqLength));
			out.write((char*)&full.matrix[0], sizeof(DNALength)* bwtSeqLength * AlphabetSize);
		}
	}

	int Read(istream &in, int _hasDebugInformation) {
		hasDebugInformation = _hasDebugInformation;
		in.read((char*) &majorBinSize, sizeof(majorBinSize));
		in.read((char*) &minorBinSize, sizeof(minorBinSize));
		in.read((char*) &numMajorBins, sizeof(numMajorBins));
		if (numMajorBins > 0) {
			major.Resize(numMajorBins * AlphabetSize);
			in.read((char*) major[0], sizeof(T_Major) * numMajorBins * AlphabetSize);
			major.nRows = numMajorBins;
			major.nCols = AlphabetSize;
		}
		in.read((char*) &numMinorBins, sizeof(numMinorBins));
		if (numMinorBins > 0) {
			minor.Resize(numMinorBins * AlphabetSize);
			in.read((char*) minor[0], sizeof(T_Minor) * numMinorBins*AlphabetSize);
			minor.nRows = numMinorBins;
			minor.nCols = AlphabetSize;
		}
		if (hasDebugInformation) {
			DNALength bwtSeqLength;
			in.read((char*)&bwtSeqLength, sizeof(bwtSeqLength));
			full.matrix = new DNALength[bwtSeqLength *AlphabetSize];
			full.nRows = bwtSeqLength;
			full.nCols = AlphabetSize;
			in.read((char*)&full.matrix[0], sizeof(DNALength)* bwtSeqLength * AlphabetSize);
		}
		return 1;
	}

};


#endif
