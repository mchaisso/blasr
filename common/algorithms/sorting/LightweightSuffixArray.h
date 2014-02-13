#ifndef ALGORITHMS_SORTING_LIGHTWEIGHT_SUFFIX_ARRAY_H_
#define ALGORITHMS_SORTING_LIGHTWEIGHT_SUFFIX_ARRAY_H_

#include <algorithm>
#include "qsufsort.h"
#include "MultikeyQuicksort.h"
#include "DifferenceCovers.h"
#include "../../Types.h"
using namespace std;

/*
 * a - b potentially may not fit into a signed type.  Use some logic
 * to get around that.
 */

UInt DiffMod(UInt a, UInt b, UInt d) {
	if (b > a) {
		return (d - ((b - a) % d))%d;
	}
	else {
		return (a - b) % d;
	}
}

void BuildDiffCoverReverseLookup(UInt diffCover[], UInt diffCoverLength,
																 UInt reverseDiffCover[] // of size diffCoverSize
																 ){
	UInt i;
	for (i = 0; i < diffCoverLength; i++ ){
		reverseDiffCover[diffCover[i]] = i;
	}
}

UInt DiffCoverFindH(UInt diffCover[], UInt diffCoverLength, UInt diffCoverSize, UInt textSize) {
	UInt h;
	for (h = 0; h < diffCoverSize; h++) {
		UInt rem = textSize % diffCoverSize ;
		if (rem == 0) return 0;
		if ((h < diffCoverSize -1 and (diffCover[h] <= rem and rem < diffCover[h+1])) or
				(h == diffCoverSize-1 and (diffCover[h] <= rem and rem < diffCoverSize))) {
			return h;
		}
	}
	return h;
}


class DiffCoverMu {
 public:
	UInt *diffCoverReverseLookup;
	UInt diffCoverLength;
	UInt diffCoverSize;
	UInt textSize;
	UInt h;
	UInt *diffCover;

	UInt compute(UInt i, UInt j) {
		return textSize/diffCoverSize * i + min(i,h+1) + j;
	}
	UInt operator()(const UInt k) {
		//
		// k is from 0 .. n (size of string)
		//
		UInt di = k % diffCoverSize;
		UInt j  = k / diffCoverSize;
		UInt i  = diffCoverReverseLookup[di];
		//		return (textSize/diffCoverSize)*i + min(i,h) + j;
		//		return (textSize/diffCoverSize)*i + i + j;
		UInt itemsInBucket; 
		//		return min(i,h)*(1 + textSize / diffCoverSize) + (i > h ? i - h : 0)*(textSize/diffCoverSize) + j;
		return (textSize/diffCoverSize) * i + min(i,h+1) + j;
	}

	DiffCoverMu() {
		diffCoverReverseLookup = NULL;
    diffCoverLength = diffCoverSize = textSize = h = 0;
    diffCoverReverseLookup = diffCover = NULL;
	}

	~DiffCoverMu() {
		if (diffCoverReverseLookup !=NULL) 
			delete[] diffCoverReverseLookup;
	}
  
	void Initialize(UInt diffCoverP[], UInt diffCoverLengthP, UInt diffCoverSizeP, UInt textSizeP) {
		diffCoverReverseLookup = new UInt[diffCoverSizeP];
		diffCoverLength = diffCoverLengthP;
		textSize        = textSizeP;
		diffCoverSize   = diffCoverSizeP;
		diffCover       = diffCoverP;
		UInt i;
		for (i = 0; i < diffCoverSize; i++ ) {
			diffCoverReverseLookup[i] = 9999999;
		}
		BuildDiffCoverReverseLookup(diffCoverP, diffCoverLength, diffCoverReverseLookup);
		h               = DiffCoverFindH(diffCoverP, diffCoverLength, diffCoverSize, textSize);
	}

};

void BuildDiffCoverLookup(UInt diffCover[], UInt diffCoverLength, UInt v, UInt diffCoverLookup[]) {
	UInt h;
	// Initialize with sentinal that shows a value has not been set (for small problems);
	for (h = 0; h < v; h++) {
		diffCoverLookup[h] = 99999999;
	}
	for (h = 0; h < v; h++) {
		//
		// For now, fill table via exhaustive search.
		//
		UInt hd;
		for (hd = 0; hd < diffCoverLength; hd++) {
			UInt dcm = (diffCover[hd] + h ) % v;
			UInt hdi;
			for (hdi = 0; hdi < diffCoverLength; hdi++ ){
				if (dcm == diffCover[hdi])
					break;
			}
			if (hdi < diffCoverLength) {
				diffCoverLookup[h] = diffCover[hd];
				break;
			}
		}
	}
}

class DiffCoverDelta {
 public:
	UInt *diffCoverLookup;
	UInt diffCoverSize;
	
	void Initialize(UInt diffCoverP[], UInt diffCoverLengthP, UInt diffCoverSizeP) {
		diffCoverLookup = new UInt[diffCoverSizeP];
		diffCoverSize   = diffCoverSizeP;
		BuildDiffCoverLookup(diffCoverP, diffCoverLengthP, diffCoverSizeP, diffCoverLookup);
	}

	UInt operator()(UInt i, UInt j) {
		return DiffMod(diffCoverLookup[DiffMod(j,i,diffCoverSize)],i,diffCoverSize);
	}
	
};



UInt NCompareSuffices(unsigned char text[], UInt a, UInt b, UInt n) {
	// not sure how to make lower amortized cost of the comparison.
	return(strncmp((const char*) &text[a], (const char*)&text[b], n));
}

UInt ComputeDSetSize(UInt diffCover, UInt diffCoverLength, UInt diffCoverSize, UInt textSize) {
	UInt div = textSize / diffCoverSize + 1;
	UInt rem = textSize % diffCoverSize;
	return div*diffCoverSize + rem;
}

void ComputeSufVNaming(UInt diffCover[], UInt diffCoverLength, UInt diffCoverN, UInt textSize, UInt lexNaming[], 
											 DiffCoverMu &mu,
											 UInt sufVNaming[]) {
	UInt nDiffCover = textSize / diffCoverN + 1;
	UInt cover;
	UInt d;
	UInt diffCoverIndex;
	UInt ln = 0;
	for (cover = 0; cover < nDiffCover; cover++) {
		for (d = 0; d < diffCoverLength; d++) {
			diffCoverIndex = cover * diffCoverN + diffCover[d];
			sufVNaming[mu(diffCoverIndex)] = lexNaming[ln];
			ln++;
		}
	}
}

UInt IndexToDiffCoverIndex(UInt index, UInt diffCoverlookup[], UInt diffCoverSize, UInt diffCoverLength ){
	UInt diff = index / diffCoverSize;
	UInt offset = index % diffCoverSize;
	return diff*diffCoverLength + diffCoverlookup[offset];
}

void DiffCoverComputeLOrder(UInt sufVNaming[], UInt sufVNamingLength, UInt maxVNaming, UInt textLength,
														DiffCoverMu &mu,
														UInt lOrder[]) {
	//
	// the sufvnaming now contains the 
	UInt i, di;
	UInt nDiffCover = textLength / mu.diffCoverSize + 1;
	UInt dci;
	for (i = 0; i < sufVNamingLength; i++ ) {
		lOrder[i] = 0;
	}

	for (dci = 0; dci < nDiffCover; dci++ ) {
		for (di = 0; di < mu.diffCoverLength; di++ ){
			i = dci*mu.diffCoverSize + mu.diffCover[di];
			if (i >= textLength) {
				break;
			}
			UInt dsetIndex = IndexToDiffCoverIndex(i, mu.diffCoverReverseLookup, mu.diffCoverSize, mu.diffCoverLength);
			UInt mui = mu(i);
			lOrder[mui] = sufVNaming[dsetIndex] + 1;
		}
	}
	lOrder[sufVNamingLength] = 0;
	//
	// The result of the sufsort function is to store the inverse suffix
	// array in the first parameter.
	//

	LarssonSuffixSort<UInt> sufsorter;
	sufsorter.INDEX_MAX = maxVNaming + 2;
	sufsorter(lOrder, sufVNaming, sufVNamingLength, maxVNaming+2, 0);
	for (i = 0; i < sufVNamingLength; i++) {
		assert(lOrder[i] > 0);
		lOrder[i]--;
	}
}


	
/*
 * Build the lex naming of the v-ordered suffices.  
 *
 * Input: textVOrder - the v-ordering of a subset of the text.
 *        textSize   - the size of the v-order set.
 *        diffCover  - the diff cover used, and it's length
 *        diffCoverLength 
 *        diffCoverSize - the size of the diff cover.
 * Output: lexNaming: the lex-naming of the v-order suffices.  The
 *        names are implemented as unsigned integers. 
 * Returns: the largest value of the lex-ordering.
 */
UInt DiffCoverBuildLexNaming( unsigned char text[], UInt textSize,
															UInt textVOrder[],
															UInt dSetSize, UInt diffCover[], UInt diffCoverLength, UInt diffCoverSize, 
															UInt diffCoverLookup[],
															UInt lexNaming[]) {
	UInt nCovers = textSize / diffCoverSize + 1;
	UInt cover = 0;
	UInt d;
	UInt vOrder = 0;
	UInt prevCoverIndex, coverIndex;
	UInt lexOrder = 0;
	UInt lexIndex = 0;
	//
	// Make sure there is something to do here.
	//
	if (dSetSize == 0) 
		return 0;
	UInt dcindex;
	dcindex = IndexToDiffCoverIndex(textVOrder[0], diffCoverLookup, diffCoverSize, diffCoverLength);
	lexNaming[dcindex] = 0;
	for (d = 1; d < dSetSize; d++) {
		if (NCompareSuffices(text, textVOrder[d-1], textVOrder[d], diffCoverSize) != 0) {
			lexOrder++;
		}
		dcindex = IndexToDiffCoverIndex(textVOrder[d], diffCoverLookup, diffCoverSize, diffCoverLength);
		lexNaming[dcindex]  = lexOrder;
	}
	return lexOrder;
}

class DiffCoverCompareSuffices {
 public:
	UInt *lOrder;
	DiffCoverDelta *delta;
	UInt diffCoverSize;
	UInt diffCoverLength;
	UInt *diffCoverReverseLookup;
	int operator()(UInt a, UInt b) {
		UInt aDCIndex, bDCIndex;
		UInt dab = (*delta)(a,b);
		aDCIndex = IndexToDiffCoverIndex(a + dab, diffCoverReverseLookup, diffCoverSize, diffCoverLength);
		bDCIndex = IndexToDiffCoverIndex(b + dab, diffCoverReverseLookup, diffCoverSize, diffCoverLength);
		return (lOrder[aDCIndex] < lOrder[bDCIndex]);
	}
};

bool LightweightSuffixSort(unsigned char text[], UInt textLength, UInt *index, int diffCoverSize) {
	//
	// index is an array of length textLength that contains all
	// suffices.
	//
	
	//
	// Phase 0. Compute delta function for difference cover.
	//
	
	//  For now, use a very small hard wired diff cover for testing
	UInt *diffCover;
	UInt diffCoverLength;
	if (InitializeDifferenceCover(diffCoverSize, diffCoverLength, diffCover) == 0) {
		cout << "ERROR! There is no difference cover of size " << diffCoverSize << " that is precomputed." << endl;
		exit(1);
	}

	DiffCoverDelta delta;

	delta.Initialize(diffCover, diffCoverLength, diffCoverSize);
	
	//
	// Phase 1. Sort suffices whose starting position modulo v is in D.
	//
	
	// The set d is given by 
	// Step 1.1 v-sort D-sample suffices
	UInt dIndex = 0; // index in D-sample
	UInt tIndex = 0; // index in text
	UInt nDiffCover;
	nDiffCover = textLength / diffCoverSize + 1;
	UInt coverIndex;
	UInt d;
	bool done = false;

	for (coverIndex = 0; coverIndex < nDiffCover and done == false; coverIndex++) {
		for (d = 0; d < diffCoverLength and done == false; d++) {
			tIndex = coverIndex * diffCoverSize + diffCover[d];
			if (tIndex >= textLength) {
				done = true;
				break;
			}
			index[dIndex] = tIndex;
			dIndex++;
		}
	}
	UInt dSetSize = dIndex;
	cerr << "Sorting " << diffCoverSize << "-prefixes of the genome." << endl;
	MediankeyBoundedQuicksort(text, index, dIndex, 0, dSetSize, 0, diffCoverSize);
	UInt i;
	
	//
	// Step 1.2 Compute l^v(i) for all i \in D_n by traversing the
	// D-sample suffixes in lexicographic order and construct s^\prime
	// by setting s^\prime[\mu(i)] = l^v(i)
	//
	UInt *lexVNaming;
	lexVNaming = new UInt[dSetSize+1];
	if (lexVNaming == NULL) {
		cout << "Could not initialize welterweight order structure." << endl;
		exit(1);
	}
	DiffCoverMu mu;
	mu.Initialize(diffCover, diffCoverLength, diffCoverSize, textLength);
	UInt largestLexName;
	cerr << "Enumerating " << diffCoverSize << "-prefixes." << endl;
	largestLexName = DiffCoverBuildLexNaming(text, textLength,
        index, dSetSize, diffCover, diffCoverLength, diffCoverSize, 
		mu.diffCoverReverseLookup, lexVNaming);
	//
	// Step 1.3 Compute ISA' of lex-order.
	//
	
	UInt *lexOrder;
	//
	// lexVNaming is allocated space.  The suffix sorting needs an
	// auxiliary array.  Since the index is not being used right now,
	// use that as the extra space.
	//
	UInt t;
	UInt dci,di;	
	for (dci = 0; dci < nDiffCover; dci++) {
		for (di = 0 ; di < diffCoverLength; di++) {
			i = dci*diffCoverSize + diffCover[di];
			if (i > textLength) {
				break;
			}
			mu.compute(di, dci);
		}
	}

	UInt *tmpLexOrder = index; 
	DiffCoverComputeLOrder(lexVNaming, dSetSize, largestLexName, textLength, mu, tmpLexOrder); 	
	lexOrder = lexVNaming;

	nDiffCover = textLength / diffCoverSize + 1;
	
	for (dci = 0; dci < nDiffCover; dci++) {
		for (di = 0 ; di < diffCoverLength; di++) {
			i = dci*diffCoverSize + diffCover[di];
			if (i >= textLength) {
				break;
			}
			UInt lexOrderIndex = mu(i);
			UInt diffCoverIndex = IndexToDiffCoverIndex(i, mu.diffCoverReverseLookup, diffCoverSize, diffCoverLength);
			lexOrder[diffCoverIndex] = tmpLexOrder[lexOrderIndex];
		}
	}

	//
	// Phase 2. Construct SA by exploiting the fact that for any i,j\in
	// [0,n-v], the relative order of the suffixes starting at
	// i+\delta(,j) and j+\delta(i,j) is already known.
	//
	cerr << "Sorting suffices." << endl;
	// Step 2.1 v-order suffices using multikey quicksort
	for (i = 0; i < textLength; i++ ){
		index[i] = i;
	}
	MediankeyBoundedQuicksort(text, index, textLength, 0, textLength, 0, diffCoverSize);

	// Step 2.2. For each group of suffixes that remains unsorted
	// (shares a prefix of length diffCoverSize, complete the sorting
	// with a comparison based on the sorting algorithm using
	// l(i+\delta(i,j)) nad l(j+\delta(i,j)) as keys when comparing
	// suffixes S_i and S_j.
	//
	DiffCoverCompareSuffices lOrderComparator;
	lOrderComparator.lOrder = lexOrder;
	lOrderComparator.delta  = &delta;
	lOrderComparator.diffCoverSize = diffCoverSize;
	lOrderComparator.diffCoverLength=diffCoverLength;
	lOrderComparator.diffCoverReverseLookup = mu.diffCoverReverseLookup;
	UInt setBegin, setEnd;
	setBegin = setEnd = 0;
	cerr << "Sorting buckets." << endl;
	int percentDone = 0;
	int curPercentage = 0;
	while(setBegin < textLength) {
		setEnd = setBegin;
		percentDone = (int)(((1.0*setBegin) / textLength) * 100);
		if ( percentDone > curPercentage) {
			cerr << " " << percentDone << "% of buckets sorted."  << endl;
			curPercentage = percentDone;
		}
		while(setEnd < textLength and
					NCompareSuffices(text, index[setBegin], index[setEnd], diffCoverSize) == 0) {
			setEnd++;
		}
		std::sort(&index[setBegin], &index[setEnd], lOrderComparator);
		setBegin = setEnd;
	}
	
    return true;
	// DONE!!!!!

}

#endif
