#ifndef RADIX_SORT_H_
#define RADIX_SORT_H_

#include <vector>
#include <iostream>
#include <assert.h>

#include "NucConversion.h"
#include "DNASequence.h"
#include "defs.h"

using namespace std;
inline
void RadixSortRL(int v[], int vs[], Nucleotide seq[], int length, int keyLength, int alphSize ) {
	// 
	// Perform an in-place radix sort on V. RL means that 
	// the least-significant digit is on the right of the key, and
	// most on left.  That's big endian, but easier to remember.
	//
	
	int *c = new int[alphSize];
	int d;

	int i;
	for (d = keyLength - 1; d >= 0; d--) {
		//
		// Count the digits.
		//
		std::fill(c, c+alphSize, 0);
		for (i = 0; i < length; i++) { c[seq[v[i] + d]]++; }
		
		//
		// Initialze the next-place buckets.
		//

		for (i = 1; i < alphSize; i++) { c[i] = c[i] + c[i-1];}
		//
		// Do value-swapping to get items in order of their bucket.
		//
		for (i = length - 1; i >= 0; i-- ) {
			//			cout << c[seq[v[i]+d]] << " gets " << v[i] << endl;
			vs[c[seq[v[i]+d]]-1] = v[i];
			c[seq[v[i]+d]]--;
		}
		// swap sorted and work arrays.
		char *key = new char[alphSize + 1];
		key[alphSize] = '\0';

		for(i = 0;i<length;i++) {
			v[i] = vs[i];
		}

	}												 
	delete[] c;
}


inline
int RadixSortLR(int v[], Nucleotide seq[], int length, int keyLength, int alphSize, int digit=0 ) {
	// 
	// Perform an in-place radix sort on V. LR for processing
	// strings from left to right (Most Significant Digit).
	// 
	//

	if (length <= 1) 
		return 0;

	if (digit >= keyLength) 
		return 1;

	int *c = new int[alphSize];
	int *b = new int[alphSize];
	int i;
	
	//
	// Count the digits.
	//
	std::fill(c, c+alphSize, 0);
	for (i = 0; i < length; i++) { c[seq[v[i] + digit]]++; }
		
	//
	// Initialze the next-place buckets.
	//
	b[0] = 0;
	for (i = 1; i < alphSize; i++) { b[i] = c[i-1] + b[i-1];}
	for (i = 1; i < alphSize; i++) { c[i] = c[i] + c[i-1];}
	
	//
	// Do value-swapping to get items in order of their bucket.
	//
	int ch = 0;
	char cur;
	i = 0;
	while (i < length) {
		cur = seq[v[i]+digit];
		if (cur == ch) {
			// cur is in the correct slot.
			i++;
		}
		else {
			// cur char is not in the correct slot, Swap it into the
			// correct slot, but don't advance i since the swap
			// may also not bring in the correct value.

			while (seq[v[b[cur]] + digit ] == cur) b[cur]++;

			// now b[cur] points to a value that does not have car 'cur' at
			// digit d, but it is in the bucket for 'cur'.
			int dest = b[cur];
			int temp = v[dest];
			v[dest] = v[i];
			v[i] = temp;
		}
		while (i == c[ch]) ch++;
	}
	RadixSortLR(&v[0], seq, c[0], keyLength, alphSize, digit+1);
	for (i = 1; i < alphSize; i++ ) {
		RadixSortLR(&v[c[i-1]], seq, c[i] - c[i-1], keyLength, alphSize, digit + 1);
	}

	delete[] b;
	delete[] c;
	return 1;
}


inline
int RadixSuffixSortLR(int v[], Nucleotide seq[], int seqLength, int keyLength, int alphSize, int digit=0, int seqEnd=0 ) {
	// 
	// Perform an in-place radix sort on V. LR for processing
	// strings from left to right (Most Significant Digit).
	//
	// By sorting on suffixes, this assumes that x$ < x*, where $ represents the end of the string, and $ 
	// is lex-less than any other character.
	// 

	//	cout << "starting out at length: " << seqLength << endl;

	if (seqLength <= 1) 
		return 0;

	if (digit >= keyLength) 
		return 1;

	if (seqEnd == 0)
		seqEnd = seqLength;

	int *c = new int[alphSize];
	int *b = new int[alphSize];
	int *boundary = new int[alphSize+1];
	int i;
	
	//
	// Count the digits.
	//
	std::fill(c, c+alphSize, 0);
	int nEnd = 0;
	int endPos = 0;
	for (i = 0; i < seqLength; i++) {
		if (v[i] + digit < seqEnd) c[seq[v[i] + digit]]++;
		else {
			//			cout << "found end pos at: " << i << endl;
			nEnd++;
			endPos = i;
		}
	}
	assert(nEnd <= 1);
		
	//
	// Initialze the next-place buckets.  b[i] implies that the place to put 
	//
	b[0] = nEnd;
	boundary[0] = nEnd;
	for (i = 0; i < alphSize; i++) { boundary[i+1] = boundary[i] + c[i];}
	for (i = 1; i < alphSize; i++) { b[i] = c[i-1] + b[i-1];}
	c[0] += nEnd;
	for (i = 1; i < alphSize; i++) { c[i] = c[i] + c[i-1];}

	if (nEnd) {
		// 
		// The suffix passes past the end of the string. 
		// Since end of string is lex-less than anything else, promote it to the
		// beginning.
		//		cout << "putting: " << v[endPos] << " at beginning of buckets." << endl;
		SWAP(v[0], v[endPos]);
	}
	//
	// Do value-swapping to get items in order of their bucket.
	//
	int ch = 0;
	char cur;
	i = nEnd;
	while (i < seqLength - 1) {
		//
		// Look to see if the current suffix-keyword contains 
		// a suffix that passes past the end of the string.  If so, 
		// the keyword  
		//
		cur = seq[v[i]+digit];
		if (cur == ch) {
			// cur is in the correct slot.
			i++;
		}
		else {
			// cur char is not in the correct slot, Swap it into the
			// correct slot, but don't advance i since the swap
			// may also not bring in the correct value.

			while (seq[v[b[cur]] + digit ] == cur) b[cur]++;

			// now b[cur] points to a value that does not have car 'cur' at
			// digit d, but it is in the bucket for 'cur'.
			int dest = b[cur];
			SWAP(v[dest], v[i]);
		}
		while (i == c[ch]) ch++;
	}
	/*
	cout << "v[0]: " << v[0] << " length: " << seqLength << endl;
	char *key = new char[keyLength + 1];
	for (i = 0; i < seqLength; i++) {
		int j;
		int wordLength = keyLength;
		if (v[i] + keyLength > seqEnd) {
			wordLength = seqEnd - v[i];
		}
		key[wordLength] = '\0';
		for (j = 0; j < wordLength; j++ ) {key[j] = TwoBitToAscii[seq[v[i] + j]];}
		cout << key << endl;
	}
*/
	//	cout << "sorting from " << boundary[0] << " ... " << boundary[1] << endl;
	RadixSuffixSortLR(&v[boundary[0]], seq, boundary[1] - boundary[0], keyLength, alphSize, digit+1, seqEnd);
	for (i = 1; i < alphSize; i++ ) {
		//		cout << "sorting from " << boundary[i] << " ... " << boundary[i+1] << endl;
		RadixSuffixSortLR(&v[boundary[i]], seq, boundary[i+1] - boundary[i], keyLength, alphSize,digit+1, seqEnd);
	}

	delete[] b;
	delete[] c;
	delete[] boundary;
	return 1;
}


#endif
