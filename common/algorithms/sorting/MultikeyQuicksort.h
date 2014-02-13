#ifndef ALGORITHMS_SORTING_MULTIKEY_QUICKSORT_H_
#define ALGORITHMS_SORTING_MULTIKEY_QUICKSORT_H_

/*
 * This is an implementation of MultiKey Quicksort, or ssort1 from
 * Bentley and Sedgewick, Fast Algorithms for Sorting and Searching
 * Strings, Proc. 8th Annual ACM-SIAM Symposium on Discrete Algorithms
 * (SODA), pages 360-369, January 1997. 
 *
 * The implementation here sorts lists of indices into an array of
 * substrings rather than the substrings themselves (for use in suffix
 * array generation).  Furthermore, it is made to be sse-alizeable to
 * speed up sorts on modern architectures. 
 */

#include <limits.h>
#include <vector>
#include <algorithm>
#include "../../NucConversion.h"
#include "../../FASTASequence.h"
#include "../../Types.h"

using namespace std;
typedef unsigned int UInt;
void UIntSwap(unsigned int &a, unsigned int &b) {
	//
	// There probably is an SSE for this.
	//
  unsigned int tmp = a;
  a = b; b = tmp;
}
	
void VecSwap( UInt i, UInt j, UInt n, UInt index[]) {
	UInt ni;
	for (ni = 0; ni < n; ni++) {
		UIntSwap(index[i], index[j]);
		i++; j++;
	}
}

unsigned char ComputeMedianValue(unsigned char text[], UInt index[], int length, 
				 UInt low, UInt high, int offset, UInt maxPossible, UInt *freq ) {
  /*
   * Compute the median value of positions at suffix+offset across all
   * suffices from [low,high).
   */

  int f;
  for (f = 0; f< maxPossible+1; f++) {
    freq[f] = 0;
  }
		
  UInt pos;
  UInt maxValue = 0;
  int val;
  //  return text[index[(high+low)/2]+offset];
  // Compute frequencies of nucleotides (+N), packed into 3-bit representation.
  for (pos = low; pos < high; pos++) {
    //		if (index[pos] + offset < length) {
    //			val = ThreeBit[text[index[pos] + offset]];
    val = text[index[pos]+offset];
    if (val > maxValue) {
      maxValue = val;
    }
    freq[val]++;
    //		}
  }
  UInt medianBound = (high - low)/2;
  int i;
  int runningTotal = 0;

  for(i =1; i<= maxValue; i++ ){
    runningTotal = freq[i] + runningTotal;
    if (runningTotal >= medianBound) {
      return i;
    }
  }
  return maxValue;
}
	
UInt FindFirstOf(unsigned char text[], UInt index[], UInt length, 
								 UInt low, UInt high, int offset, Nucleotide nuc) {
	UInt p;
	for (p = low; p < high and text[index[p]+offset] != nuc; p++);
	return p;
}

inline
void SwapIndices(UInt &a, UInt &b) {
	UInt temp = a;
	a = b;
	b = temp;
}

/*
UInt min(UInt a, UInt b) {
	return (a < b) ? a : b;
}
*/
void TransformSequenceForSorting(unsigned char text[], UInt textLength, int bound) {
	UInt i;
	for (i = 0; i < textLength; i++) {
		text[i] = ThreeBit[text[i]] + 1;
	}
	for (i = textLength; i < textLength+bound;i++) {
		text[i] = 0;
	}
}

void TransformBackSequence(Nucleotide text[], UInt textLength) {
	UInt i;
	for (i = 0; i < textLength; i++) {
		text[i]--;
		text[i] = ThreeBitToAscii[text[i]];
	}
}	


/*
 * depth: the current depth of how much is sorted.
 * bound: how far to sort.
 */
void MediankeyBoundedQuicksort(unsigned char text[], UInt index[], UInt length,
			       UInt low, UInt high, int depth, int bound, UInt maxChar= 0, UInt *freq=NULL) {
	
	if (high - low <= 1) 
		return;
	if (depth > bound)
		return;

	if (freq == NULL) {
	  UInt ci;
	  maxChar = 0;
	  for (ci = low; ci < high; ci++) {
	    UInt c = text[index[ci]+depth];
	    if ( c > maxChar) {
	      maxChar = c;
	    }
	  }
	  freq = new UInt[maxChar+1];
	}

	Nucleotide medianChar = ComputeMedianValue(text, index, length, low, high, depth, maxChar, freq );
	UInt medianCharPos    = FindFirstOf(text, index, length, low, high, depth, medianChar);
	UInt medLeft, lastLeft;
	UInt medRight, lastRight;
	
	// 
	// Pack the median characters into the sides of the array.
	//
	SwapIndices(index[low], index[medianCharPos]);
	
	medLeft = lastLeft = low+1;
	medRight = lastRight = high-1;

	for(;;) {
		Nucleotide nuc;
		while(lastLeft <= lastRight and (nuc = text[index[lastLeft]+depth]) <= medianChar) {
			if (nuc == medianChar) {
				SwapIndices(index[medLeft], index[lastLeft]);
				medLeft++;
			}
			lastLeft++;
		}

		while(lastLeft <= lastRight and (nuc = text[index[lastRight]+depth]) >= medianChar) {
			if (nuc == medianChar) {
				SwapIndices(index[medRight], index[lastRight]);
				medRight--;
			}
			lastRight--;
		}
		if (lastLeft > lastRight) {
			// done with median sort.
			break;
		}
		//
		// Otherwise, this ends with an index at the left out of order
		// from the right
		assert(text[index[lastLeft]+depth] > text[index[lastRight]+depth]);
		SwapIndices(index[lastLeft], index[lastRight]);
	}
	
	//
	// Now join the indices of the median charactes in the middle of the
	// array. 
	//
	// move left outside to middle
	UInt swapLeftLength = min(medLeft - low, lastLeft - medLeft);
	VecSwap(low, lastLeft - swapLeftLength, swapLeftLength, index);
	// move right outside to middle
	UInt swapRightLength = min(high - medRight - 1, medRight - lastRight);
	VecSwap(lastRight+1, high-swapRightLength, swapRightLength, index);
	
	UInt medianStartIndex = low + lastLeft - medLeft;
	UInt medianEndIndex   = lastRight + (high - medRight);
	
	// Sort the suffices with keys lower than the median. 
	// Since these may contain multiple keys that are less than the
	// median, the same depth is used.
	MediankeyBoundedQuicksort(text, index, length, 
				  low, 
				  medianStartIndex,
				  depth, bound, maxChar,freq);
	if (medianEndIndex - medianStartIndex > 1) {
	  MediankeyBoundedQuicksort(text, index, length,
				    medianStartIndex, medianEndIndex, depth+1, bound, maxChar, freq);
	}
	
	MediankeyBoundedQuicksort(text, index, length,
				  medianEndIndex,
				  high,
				  depth, bound, maxChar, freq);
	

}

#endif
