#ifndef SUFFIX_ARRAY_H_
#define SUFFIX_ARRAY_H_
#include <string.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "LCPTable.h"
#include "defs.h"
#include "tuples/DNATuple.h"
#include "tuples/CompressedDNATuple.h"
#include "algorithms/compare/CompareStrings.h"
#include "algorithms/sorting/qsufsort.h"
#include "algorithms/sorting/LightweightSuffixArray.h"
#include "qvs/QualityValue.h"
#include "DNASequence.h"
#include "NucConversion.h"
#include "utils/ProtectedNew.h"
/*
 * Suffix array implementation, with a Manber and Meyers sort, but
 * that is typically not used.
 *
 */


using namespace std;

typedef enum E_SAType {manmy, slow, mcilroy, larsson, kark, mafe, welter} SAType;

template<typename T>
class CompareSuffixes {
 public:
	T t;
	int refLength;
	CompareSuffixes(T tref, int prefLength) {
		t = tref;
		refLength = prefLength;
	}
	int operator()(int a, int b) {
		int aSufLen = refLength - a;
		int bSufLen = refLength - b;
		int abMinLength = MIN(aSufLen, bSufLen);
		int cmpRes = memcmp(&(t[a]), &(t[b]), abMinLength);
		if (cmpRes == 0) {
			if (aSufLen < bSufLen) {
				return 1;
			}
			else {
				return 0;
			}
		}
		else {
			return cmpRes < 0;
		}
	}
};

typedef uint32_t SAIndex;
typedef uint32_t SAIndexLength;

template<typename T, 
	typename Sigma,
	typename Compare = DefaultCompareStrings<T>,
	typename Tuple   = DNATuple >
	class SuffixArray {
 public:
 SAIndex *index;
 bool deleteStructures;
 T*  target;
 SAIndex length;
 SAIndex *startPosTable, *endPosTable;
 SAIndexLength lookupTableLength;
 SAIndex lookupPrefixLength;
 TupleMetrics tm;
 unsigned int magicNumber;
 unsigned int ckMagicNumber;
 typedef Compare CompareType;
 enum Component { CompArray, CompLookupTable, CompLCPTable};
 static const int ComponentListLength = 2;
 static const int FullSearch = -1;

 SAIndex operator[](SAIndex i) {
	 return index[i];
 }
 int componentList[ComponentListLength];

 // vector<SAIndex> leftBound, rightBound;

 inline	int LengthLongestCommonPrefix(T *a, int alen, T *b, int blen) {
	 int i;
	 for (i = 0 ; i < alen and i < blen; i++ ) 
		 if (a[i] != b[i])
			 break;
	 return i;
 }

 SuffixArray() {
	 // Not necessarily using the lookup table.
   // The magic number is linked with a version 
	 magicNumber = 0xacac0001;
	 startPosTable = endPosTable = NULL;
	 lookupPrefixLength = 0;
	 lookupTableLength = 0;
	 deleteStructures  = true;
     ckMagicNumber = 0;
   length = 0;
   int i;
   for (i = 0; i < ComponentListLength; i++) {
     componentList[i] = false;
   }
	 // Must create a suffix array, but for now make it null.
	 target = NULL;
     index = NULL;
 }
 ~SuffixArray() {
	 if (deleteStructures == false) {
		 //
		 // It is possible this class is referencing another structrue. In
		 // this case, do not try and delete in the destructor.
		 //
		 return;
	 }
	 if (startPosTable != NULL) {
		 delete[] startPosTable;
	 }
	 if (endPosTable != NULL) {
		 delete[] endPosTable;
	 }
	 if (index != NULL) {
		 delete[] index;
	 }
 }

 int StringLessThanEqual(T *a, int aLen, T *b, int bLen) {
	 return Compare::LessThanEqual(a, aLen, b, bLen);
 }

 int StringEquals(T *a, int aLen, T *b, int bLen) {
	 return Compare::Equal(a, aLen, b, bLen);
 }
 int StringLessThan(T *a, int aLen, T *b, int bLen) {
	 return Compare::LessThan(a, aLen, b, bLen);
 }
	
 void InitAsciiCharDNAAlphabet(vector<int> &dnaAlphabet) {
	 int i;
	 for (i = 0; i < 127; i++) {
		 dnaAlphabet.push_back(i);
	 }
 }

 void InitTwoBitDNAAlphabet(vector<int> &dnaAlphabet) {
	 dnaAlphabet.push_back(0);
	 dnaAlphabet.push_back(1);	
	 dnaAlphabet.push_back(2);
	 dnaAlphabet.push_back(3);
 }

 void InitThreeBitDNAAlphabet(vector<int> &dnaAlphabet) {
	 //
	 // This is initialized to have ACTG-0123, N=4, and EOF=5
	 //
	 dnaAlphabet.push_back(0);
	 dnaAlphabet.push_back(1);	
	 dnaAlphabet.push_back(2);
	 dnaAlphabet.push_back(3);
	 dnaAlphabet.push_back(4);
	 dnaAlphabet.push_back(5);
 }
	
 void PrintSuffices(T *target, int targetLength, int maxPrintLength) {
	 string seq;
	 seq.resize(maxPrintLength+1);
	 SAIndex i, s;
	 seq[maxPrintLength] = '\0';
	 for (i = 0; i < length; i++) {
		 DNALength suffixLength = maxPrintLength;
		 if (index[i] + maxPrintLength > length) {
			 suffixLength = length - index[i];
		 }
		 cout << index[i] << " " << suffixLength << " ";
		 seq.resize(suffixLength);
		 for (s = 0; s < suffixLength; s++ ){ 
			 seq[s] = TwoBitToAscii[target[index[i] + s]];
		 }
		 seq[suffixLength] = '\0';
		 cout << seq << endl;
	 }
 }

 void BuildLookupTable(T *target, SAIndexLength targetLength, int prefixLengthP) { 
		
	 //
	 // pprefixLength is the length used to lookup the index boundaries
	 // given a string.
	 //
		
	 SAIndexLength i;
	 tm.tupleSize = lookupPrefixLength = prefixLengthP;
	 tm.InitializeMask();
	 lookupTableLength = 1 << (2*lookupPrefixLength);
	 startPosTable = new SAIndex[lookupTableLength];
	 endPosTable   = new SAIndex[lookupTableLength];
	 Tuple curPrefix, nextPrefix;
	 SAIndex tablePrefixIndex = 0;
		
	 for (i = 0; i < lookupTableLength; i++) {
		 startPosTable[i] = endPosTable[i] = 0;
	 }
	 i = 0;
	 VectorIndex tablePos;
	 SAIndex     indexPos;
	 indexPos = 0;
	 do {
		 // Advance to the first position that may be translated into a tuple.
		 if (targetLength < lookupPrefixLength)
			 break;
		 while(indexPos < targetLength - lookupPrefixLength + 1 and 
					 index[indexPos] + lookupPrefixLength > targetLength) {
			 indexPos++;
		 }
		 if (indexPos >= targetLength - lookupPrefixLength + 1) {
			 break;
		 }
		 while (indexPos < targetLength - lookupPrefixLength + 1 and
						curPrefix.FromStringLR((Nucleotide*) &target[index[indexPos]], tm) == 0) {
			 ++indexPos;
		 }
		 
		 startPosTable[curPrefix.tuple] = indexPos;
		 indexPos++;
		 while(indexPos < targetLength - lookupPrefixLength + 1 and
					 index[indexPos] + lookupPrefixLength < targetLength) {
			 nextPrefix.tuple = 0;
			 nextPrefix.FromStringLR((Nucleotide*) &target[index[indexPos]], tm);
			 if (nextPrefix.tuple != curPrefix.tuple) {
				 break;
			 }
			 else {
				 indexPos++;
			 }
		 }
		 endPosTable[curPrefix.tuple] = indexPos;
	 }
	 while (indexPos < targetLength - lookupPrefixLength + 1 and 
					curPrefix.tuple < lookupTableLength-1);
 }
 
 void AllocateSuffixArray(SAIndexLength stringLength) {
	 index = ProtectedNew<SAIndex>(stringLength + 1);
	 length = stringLength;
 }

 void LarssonBuildSuffixArray(T* target, SAIndexLength targetLength, Sigma &alphabet) {
	 index =  ProtectedNew<SAIndex>(targetLength+1);
	 SAIndex *p = ProtectedNew<SAIndex>(targetLength+1);
	 SAIndexLength i;
	 for (i = 0; i < targetLength; i++) { index[i] = target[i] + 1;}
	 SAIndexLength maxVal = 0;
	 for (i = 0; i < targetLength; i++) { maxVal = index[i] > maxVal ?  index[i] : maxVal;}
	 index[targetLength] = 0;
	 LarssonSuffixSort<SAIndex, UINT_MAX> sorter;
	 sorter(index, p, ((SAIndex) targetLength), ((SAIndex) maxVal+1), (SAIndex) 1 );
	 for (i = 0; i < targetLength; i++ ){ index[i] = p[i+1];};
	 length = targetLength;
	 delete[] p;
 }

 void LightweightBuildSuffixArray(T*target, SAIndexLength targetLength, int diffCoverSize=2281) {
	 index = ProtectedNew<SAIndex>(targetLength+1);
	 length = targetLength;
	 DNALength pos;
	 for (pos = 0; pos < targetLength; pos++) {
		 target[pos]++;
	 }
	 LightweightSuffixSort(target, targetLength, index, diffCoverSize);
	 for (pos = 0; pos < targetLength; pos++) {
		 target[pos]--;
	 }
	 
 }

 void MMBuildSuffixArray(T* target, SAIndexLength targetLength, Sigma &alphabet) {
	 /*
		* Manber and Myers suffix array construction.
		*/
	 length = targetLength;
	 VectorIndex a;
	 vector<int> prm;
	 vector<int> bucket;
	 vector<int> count;
	 // To be changed to bit vectors
	 vector<bool> bh, b2h;
	 bucket.resize(alphabet.size());
		
	 prm.resize(targetLength);
	 count.resize(targetLength);
	 bh.resize(targetLength+1);
	 b2h.resize(targetLength+1);
	 std::fill(bh.begin(), bh.end(), false);
	 std::fill(b2h.begin(), b2h.end(), false);
	 std::fill(count.begin(), count.end(), 0);
	 index = new SAIndex[targetLength];
	 for (a = 0; a < alphabet.size(); a++ ) {
		 bucket[a] = -1;
	 }

	 SAIndexLength i;
	 for (i = 0; i < targetLength; i++) {
		 index[i] = bucket[target[i]];
		 bucket[target[i]] = i;
	 }
		 
	 int j;
	 SAIndex c;
	 std::fill(prm.begin(), prm.end(), -1);
	 //
	 // Prepare the buckets.
	 //
	 c = 0;
	 int b;
	 for (a = 0; a < alphabet.size(); a++) { 
		 b = bucket[alphabet[a]]; // position of last suffix starting with 'a'
		 while (b != -1) {
			 j = index[b];
			 prm[b] = c;
			 if (b == bucket[a]) {
				 bh[c] = true;
			 }
			 else {
				 bh[c] = false;
			 }
			 c = c + 1;
			 b = j;
		 }
	 }
	 b2h[targetLength] = bh[targetLength] = true;
	 // fill the index with positions sorted by the first character.
	 for (i = 0; i < targetLength; i++) {
		 index[prm[i]] = i;
	 }

	 SAIndex h;
	 h = 1;
	 SAIndex l, r;

	 while (h < targetLength) {
		 // re-order the buckets;
		 l = 0;
		 int bstart;
		 while (l < targetLength) {
			 bstart = l;
			 r = l + 1;
			 count[l] = 0;
			 //				bh[l] = 0;
			 while (bh[r] == false) {r++;} // find the begining of the next bucket.
			 while (l < r) {
				 assert(l < targetLength);
				 prm[index[l]] = bstart;
				 l++;
			 }
		 }
			
		 SAIndex d = targetLength - h;
		 SAIndex e = prm[d]; 

		 /*
			* Phase 1: Set up the buckets in the index and bh list.
			*/

		 //
		 // suffix d needs to be moved to the front of it's bucket.
		 // d should exist in the bucket starting at prm[d]
		 SAIndex i;
	
		 l = 0;
		 r = 1;

		 //
		 // Move each d that is h backwards up in it's 2h bucket.
		 //

		 d = targetLength - h; 
		 e = prm[d]; 
		 bh[e]    = true;     // e is bstart, the beginning of the bucket.
		 prm[d]   = e + count[e];
		 count[e] = count[e] + 1;
		 b2h[prm[d]] = true;

		 for (c = 0; c < targetLength; c++ ){
			 // d is T_i
			 d = index[c] - h;
			 if (index[c] >= h and d < targetLength) {
				 e           = prm[d];
				 prm[d]      = e + count[e];
				 count[e]    = count[e] + 1;
				 b2h[prm[d]] = true;
			 }
		 }


		 //
		 // Fix the bucket boundaries.
		 //

		 l = 0;


		 while(l < targetLength) {
				
			 // First assign b2h to be 1 on the entire portion of the 
			 // current bucket (from l ... bh[c]==true).
			 for (c = l; c == l or bh[c] == false; c++)  {
				 d = index[c] - h;
				 if (d >= 0 and d < targetLength) {
					 b2h[prm[d]] = true;
				 }
			 }

			 //
			 // Mark the start boundaries of the 2h bucket.
			 //
			 for (c = l; c == l or bh[c] == false; c++) {
				 d = index[c] - h;
				 if (d >= 0 and d < targetLength) {
					 if (b2h[prm[d]] == true) {
						 j = prm[d] + 1;
						 // advance j to the next bucket.
						 while (!(bh[j] == true or b2h[j] == false)) {
							 j++;
						 }

						 e = j;
						 SAIndex f;
						 for (f = prm[d] + 1; f <= e - 1; f++) { 
							 b2h[f] = false;
						 }
					 }
				 }
			 }
			 l = c;
		 }

		 for (i = 0; i < targetLength; i++) { 
			 index[prm[i]] = i;
		 }

		 for (i = 0 ; i < targetLength; i++) {
			 if (b2h[i] == true and bh[i] == false) {
				 bh[i] = b2h[i];
			 }
		 }
		 h <<= 1;
	 }
 }

 void BuildSuffixArray(T* target, SAIndex targetLength, Sigma &alphabet) {
	 length = targetLength;
	 index  = ProtectedNew<SAIndex>(length);
	 CompareSuffixes<T*> cmp(target, length);
	 SAIndex i;
	 for (i = 0; i < length; i++ ){ 
		 index[i] = i;
	 }
	 std::sort(index, index + length, cmp);
 }

 void WriteArray(ofstream &out) {
	 out.write((char*) &length, sizeof(int));
	 out.write((char*) index, sizeof(int) * (length));
 }
	
 void WriteLookupTable(ofstream &out) {

	 out.write((char*) &lookupTableLength, sizeof(SAIndex));
	 out.write((char*) &lookupPrefixLength, sizeof(SAIndex));
	 out.write((char*) startPosTable, sizeof(SAIndex) * (lookupTableLength));
	 out.write((char*) endPosTable, sizeof(SAIndex) * (lookupTableLength));
 }

void WriteComponentList(ofstream &out) {
	 //
	 // First build the component list.
	 //
	 if (index != NULL)
		 componentList[CompArray] = 1;
	 else 
		 componentList[CompArray] = 0;
		
	 if (startPosTable != NULL)
		 componentList[CompLookupTable] = 1;
	 else
		 componentList[CompLookupTable] = 0;

	 out.write((char*) componentList, sizeof(int) * ComponentListLength);
 }

 void WriteLCPTable(ofstream &out) {
	 cout << "NOT YET IMPLEMENTED." << endl;
	 exit(1);
 }

 void Write(string &outFileName) {

	 //
	 // The suffix array is written in 2 or more parts:
	 //   1 - a preamble listing the components of the
	 //       array that are written
	 //   2 - The components.
	 //
	 // 
	 ofstream suffixArrayOut;
	 suffixArrayOut.open(outFileName.c_str(), ios::binary);
	 if (!suffixArrayOut.good()) {
		 cout << "Could not open " << outFileName << endl;
		 exit(1);
	 }
	 WriteMagicNumber(suffixArrayOut);
	 // write the preamble
	 WriteComponentList(suffixArrayOut);

	 // write the components
	 if (componentList[CompArray]) {
		 WriteArray(suffixArrayOut);
	 }
	 if (componentList[CompLookupTable]) {
		 WriteLookupTable(suffixArrayOut);
	 }
	 suffixArrayOut.close();
 }
 void WriteMagicNumber(ofstream &out) {
	 out.write((char*) &magicNumber, sizeof(int));
 }

 int ReadMagicNumber(ifstream &in) {
	 in.read((char*) &ckMagicNumber, sizeof(int));
	 if (ckMagicNumber != magicNumber) {
		 return 0;
	 }
	 else { 
		 return 1;
	 }
 }

 void ReadComponentList(ifstream &in) { 
	 in.read((char*) componentList, sizeof(int) * ComponentListLength);
 }

 void ReadAllocatedArray(ifstream &in) {
	 in.read((char*) index, sizeof(int) * length);
 }

 void LightReadArray(ifstream &in) {
   in.read((char*) &length, sizeof(int));
   // skip the actual array
   in.seekg(length*sizeof(int), std::ios_base::cur);
 }

 void ReadArray(ifstream &in) {
	 in.read((char*) &length, sizeof(int));
	 index = ProtectedNew<SAIndex>(length);
	 ReadAllocatedArray(in);
 }

 void ReadAllocatedLookupTable(ifstream &in) {
	 in.read((char*) startPosTable, sizeof(int) * (lookupTableLength));
	 in.read((char*) endPosTable, sizeof(int) * (lookupTableLength));
 }

 void ReadLookupTableLengths(ifstream &in) {
	 in.read((char*) &lookupTableLength, sizeof(int));
	 in.read((char*) &lookupPrefixLength, sizeof(int));
 }

 void ReadLookupTable(ifstream &in) {
	 ReadLookupTableLengths(in);
	 tm.Initialize(lookupPrefixLength);
	 startPosTable = ProtectedNew<SAIndex>(lookupTableLength);
	 endPosTable   = ProtectedNew<SAIndex>(lookupTableLength);
	 ReadAllocatedLookupTable(in);
 }

 void ReadLCPTable(ifstream &in) {
	 cout <<" NOT YET IMPLEMENTED!!!" << endl;
	 exit(1);
 }

 bool LightRead(string &inFileName) {
	 ifstream saIn;
	 saIn.open(inFileName.c_str(), ios::binary);
	 int hasMagicNumber;
	 hasMagicNumber = ReadMagicNumber(saIn);
	 if (hasMagicNumber == 1) {
		 ReadComponentList(saIn);
     LightReadArray(saIn);
     ReadLookupTable(saIn);
     saIn.close();
     return true;
   }
   else {
     saIn.close();
     return false;
   }
 }

 bool Read(string &inFileName) {
	 ifstream saIn;
	 saIn.open(inFileName.c_str(), ios::binary);
	 int hasMagicNumber;
	 hasMagicNumber = ReadMagicNumber(saIn);
	 if (hasMagicNumber == 1) {
		 ReadComponentList(saIn);
		 if (componentList[CompArray]) {
			 ReadArray(saIn);
		 }
		 if (componentList[CompLookupTable]) {
			 ReadLookupTable(saIn);
		 }
	 saIn.close();
     return true;
	 }
	 else {
	 saIn.close();
     return false;
	 }
 }

 int SearchLCP(T* target, T* query, DNALength queryLength, SAIndex &low, SAIndex &high, DNALength &lcpLength, DNALength maxlcp) {
	 //		cout << "searching lcp with query of length: " << queryLength << endl;
	 lcpLength = 0;
	 if (startPosTable != NULL and
			 queryLength >= lookupPrefixLength) {
		 Tuple lookupTuple;
		 int left, right;
		 // just in case this was changed.
		 lookupTuple.FromStringLR(query, tm);
		 left  = startPosTable[lookupTuple.tuple];
		 right = endPosTable[lookupTuple.tuple];
		 //
		 // When left == right, the k-mer in the read did not exist in the
		 // genome.  Don't even try and map it in this case.
		 //
		 if (left == right) {
			 low = high = 0;
			 return 0;
		 }
		 
		 //
		 // Otherwise, the sequence of length 'lookupPrefixLength' was found
		 // in the genome.  The bounds of this prefix in the suffix array
		 // are stored in the lookup tables, so begin the binary search there.
		 //
		 lcpLength = lookupPrefixLength;
		 low = left, high = right;
	 }
	 else {
		 low = 0; high = length - 1;
		 lcpLength = 0;
	 }		
	 int prevLow = low;
	 int prevHigh = high;
	 int prevLCPLength = lcpLength - 1;

	 // When the boundaries and the string share a prefix, it is not necessary
	 // to use this as a comparison in further lcp searches.
	 prevLCPLength = lcpLength;
		
	 Search(target, query, queryLength, low, high, low, high, 0);

	 DNALength lowLCP = lookupPrefixLength, highLCP = lookupPrefixLength;
	 while (lowLCP < queryLength and index[low]+lowLCP < length and 
					target[index[low] + lowLCP] == query[lowLCP]) lowLCP++;

	 while (highLCP < queryLength and index[high]+highLCP < length and 
					target[index[high] + highLCP] == query[highLCP]) highLCP++;

	 DNALength minLCP = highLCP;
	 if (lowLCP < highLCP ) { 
		 minLCP = lowLCP;
	 }

	 while (minLCP >= (lookupPrefixLength -2 )and 
					low > 0 and high < (length - minLCP) and high - low < 10) {
		 while(low  > 0 and StringEquals(&target[index[low]], minLCP, &target[index[high]], minLCP)) low--;
		 while(high > 0 and StringEquals(&target[index[low]], minLCP, &target[index[high]], minLCP)) high++;
		 --minLCP;
	 }

	 //
	 // The LCP is not an exact match to the end of the string.
	 //

	 prevLow  = low;
	 prevHigh = high;

	 low = prevLow; high = prevHigh;
	 if (low < high and high - low < 100) {
		 return queryLength;
	 }
	 else {
		 high = low - 1;
		 lcpLength = 0;
	 }
	 return lcpLength;
 }

 int Search(T* target, T* query, DNALength queryLength, SAIndex left, SAIndex right, SAIndex &low, SAIndex &high, unsigned int offset=0) {
	 if (offset >= queryLength) {
		 return high - low;
	 }
	 SearchLow(target, query, queryLength, left, right, low, offset);
	 SearchHigh(target, query, queryLength, left, right, high, offset);
	 return high - low;
 }

 int Search(T* target, T* query, DNALength queryLength, SAIndex &low, SAIndex &high, int offset = 0) {

	 int left = 0;
	 int right = length - 1;
	 //
	 // Constrain the lookup if a lookup table exists.
	 //
	 if (startPosTable != NULL and
			 queryLength >= lookupPrefixLength) {
		 Tuple lookupTuple;
		 lookupTuple.FromStringLR(query, tm);
		 left  = startPosTable[lookupTuple.tuple];
		 right = endPosTable[lookupTuple.tuple];
	 }
	 return Search(target, query, queryLength, left, right, low, high, offset);
 }


 long SearchLeftBound(T* target, long targetLength, DNALength targetOffset,  T queryChar, long l, long r) {
	 long ll, lr;
	 ll = l;
	 lr = r;
	 long m;
	 long targetSufLen = 0;
	 while (ll < lr) {
		 m = (ll + lr) / 2;
		 targetSufLen = targetLength - index[m];
		 if (targetSufLen == targetOffset) {
			 ll =m + 1;
			 continue;
		 }
		 //
		 // The suffix at index[m] is shorter than the lengths of the 
		 // two sequences being compared.  With the Larsson
		 // implementation, that means that the target suffix is lex-less
		 // than the read.
		 int comp;
		 if (targetSufLen < targetOffset) {
			 comp = -1;
		 }
		 else {
			 //
			 // There is enough sequence to compare the target suffix with
			 // the query suffix.
			 //
			 assert(index[m]+targetOffset < targetLength);

			 comp = Compare::Compare(target[index[m]+targetOffset], queryChar);
       
		 }
     if (comp < 0) {
			 ll = m + 1;
		 }
		 else {
			 lr = m;
		 }
	 }
	 return ll;
 }

 long SearchRightBound(T* target, long targetLength, DNALength targetOffset, 
                       T queryChar, long l, long r) {
	 long rl, rr;
	 rl = l;
	 rr = r;
	 long m;
	 long targetSufLen;
	 while (rl < rr) {
		 m = (rl + rr) / 2;
		 targetSufLen = targetLength - index[m];
		 if (targetSufLen == targetOffset) { 
       rr = m;
			 break; 
		 }
		 if (targetSufLen < targetOffset) {
			 rr = m ;
		 }
		 else {
			 /*
				* Do not try and map stretches of N. These do not add
				* infomrative anchors.
				*/
       /*
			 if (ThreeBit[target[index[m]+targetOffset]] >= 4 or
           ThreeBit[queryChar] >= 4) {
             rl = rr;
             break;
       }
       */
			 int comp = Compare::Compare(target[index[m] + targetOffset], queryChar);
			 if (comp <= 0) {
				 rl = m + 1;
			 }
			 else {
				 rr = m ;
			 }
		 }
	 }
	 return rr;
 }

 /*
	* Search the suffix array for the bounds l and r that specify the
	* indices in the suffix array that have the longest common prefix
	* between the read and the genome.
	*/

 int SearchLCPBounds(T*target, long targetLength, T*query, DNALength queryLength, SAIndex &l, SAIndex &r, DNALength &refOffset, DNALength &queryOffset) {
	 //	 l = 0; r = targetLength;
	 for (; refOffset < targetLength and  queryOffset < queryLength and l < r; queryOffset++, refOffset++) {
		 cout << "bounds: " << l << ", " << r << endl;
		 //
		 // Band l by the character at query[offset]
		 //

		 l = SearchLeftBound(target, targetLength, refOffset, query[queryOffset], l, r);

		 //
		 // If the current search is past the end of the suffix array, it
		 // will be impossible to extend.
		 //
		 if (index[l] + refOffset >= targetLength or 
				 Compare::Compare(target[index[l] + refOffset], query[queryOffset]) != 0) {
			 break;
		 }

		 r = SearchRightBound(target, targetLength, refOffset, query[queryOffset], l, r);
		 if (Compare::Compare(query[queryOffset], target[index[l]+refOffset]) != 0 or
				 Compare::Compare(query[queryOffset], target[index[r]+refOffset]) != 0) {
			 break;
		 }
	 }
	 return refOffset;
 }


 int StoreLCPBounds(T *target, long targetLength,
										T *query,  long queryLength,
										SAIndex &low, SAIndex &high) {

	 DNALength targetOffset = 0;
	 DNALength queryOffset  = 0;
	 
	 DNALength lcpLength = 0;
	 low = 0; high = targetLength;
	 for (; index[low] + targetOffset < targetLength and
          targetOffset < targetLength  and 
					queryOffset < queryLength and 
					low < high ;
				targetOffset++, queryOffset++, lcpLength++) {
		 //
		 // Band l by the character at query[offset]
		 //

		 low = SearchLeftBound(target, targetLength, targetOffset, query[queryOffset], low, high);

		 //
		 // If the current search is past the end of the suffix array, it
		 // will be impossible to extend.
		 //
		 if (index[low] + targetOffset > targetLength or 
				 Compare::Compare(target[index[low] + targetOffset], query[queryOffset]) != 0 or
				 ThreeBit[query[queryOffset]] > 3) {
			 break;
		 }

		 high = SearchRightBound(target, targetLength, targetOffset, query[queryOffset], low, high);

	 }
	 return lcpLength;

 }

 int CountNumBranches(T* target, DNALength targetLength, DNALength targetOffset, SAIndex low, SAIndex high) {
   //
   // look to see how many different characters start suffices between
   // low and high at targetOffset
   //

   // Check some easy boundary conditions.
   //

   // 1. No branches (indices do not define any subset of the suffix
   // array). 
   if (high <= low) {
     return 0;
   }
   // 2. One branch, 
   if (target[index[low] + targetOffset ] == target[index[high-1] + targetOffset]) {
     return 1;
   }
   int numBranches = 1;
   // More than one branch.
   while ( low < high ) {
     //
     // Find the band where the suffices share the same chatacter
     // 'targetOffset' bases into the suffix as the first suffix in
     // the band given to this function.
     //
     SAIndex curCharHigh = high;
     curCharHigh = SearchRightBound(target, targetLength, targetOffset, target[index[low]+targetOffset], low, high);
     if (curCharHigh != high) {
       ++numBranches;
     }
     low = curCharHigh;
   }
   return numBranches;
 }


 int StoreLCPBounds(T *target, long targetLength, // The string which the suffix array is built on.
										T *query, DNALength queryLength, // The query string. search starts at pos 0 in this string
										bool useLookupTable,  // Should the indices of the first k bases be determined by a lookup table?
										int  maxMatchLength,  // Stop extending match at lcp length = maxMatchLength,
										// Vectors containing lcpLeft and lcpRight from 0 ... lcpLength.
										vector<SAIndex> &lcpLeftBounds, vector<SAIndex> &lcpRightBounds,
										bool stopOnceUnique=false) {

	 //
	 // Precondition: target[l][0] >= query[offset]
	 //
	 long l, r;
	 
	 l = 0; r = targetLength;
	 DNALength lcpLength = 0;
	 Tuple lookupTuple;
	 lookupTuple.tuple = -1;

	 /*
		* Various parameters may make the search through the SA not use
		* the full binary search. If priorLCP is > 0, the search for an
		* LCP is limited to lcpLeftBounds[priorLCP] and lcpRightBounds[priorLCP].
		* This is the case when continuing a search using branched
		* re-uses previous lcp searches.
		*/
	 
	 if (useLookupTable and 
			 startPosTable != NULL) {
		 // just in case this was changed.
		 if (lookupTuple.FromStringLR(query, tm)) {
			 l  = startPosTable[lookupTuple.tuple];
			 r  = endPosTable[lookupTuple.tuple];
			 lcpLength = lookupPrefixLength;
		 }
		 else {
       //
       // Not able to find a match for this sequence, so do not
       // register a hit.
       //
			 l = 0;
			 r = 0;
			 lcpLength = 0;
       return 0;
		 }
     //
     // the values of startPosTable and endPosTable are the same when
     // there are no matches.  When they are not equal, a valid range
     // has been found, so store this.
     //
		 if (l < r) {
			 VectorIndex off_i;
			 VectorIndex boundLength = lcpLeftBounds.size();
			 lcpLeftBounds.push_back(l);
			 lcpRightBounds.push_back(r);
		 }
		 else {
       //
       // No exact match found in the lookup table, do not bother
       // searching, and return 0 lcp length.
       //
			 return 0;
		 }
	 }
	 
   //
   // Search the suffix array for the longest common prefix between
   // the read and the genome.
   //
	 while( l < r and  
          lcpLength < queryLength // stop searching when the end of
                                  // the query is reached.
          ) {
     
     //
     // If there is only one match in the suffix array, and and not
     // extending matches as far as possible (stopping the search once
     // they are unique), halt the search.
     if (stopOnceUnique and l == r - 1) {
       break;
     }

     //
     // If there is a maximal match length and it is reached, stop
     // searchign as well.
     //
     if (maxMatchLength and lcpLength >= maxMatchLength) {
       break;
     }
     
     //
     // If the match extends into one or more N's, stop.  The reason
     // for this is that sometimes people will set up genome databases
     // by appending a stretch of N's between matches (although they
     // should just use a multi-fasta file).  Since the reads also
     // have stretches of N's, this tends to slow the search down
     // dramatically. 
     if (ThreeBit[target[index[l] + lcpLength]] >= 4) {
       break;
     }

		 //
		 // Find the bounds in the suffix array matching query[0... lcp]
		 // and target.
		 //

		 l = SearchLeftBound(target, targetLength, lcpLength, query[lcpLength], l, r);
		 r = SearchRightBound(target, targetLength, lcpLength, query[lcpLength], l, r);


		 //
		 // If the current search is past the end of the suffix array, it
		 // will be impossible to extend.
		 //
		 if (l == r or // if this point is reached, stop loop now since
                   // otherwise the lcp length will be incremented by
                   // 1, which will give one longer than the actual
                   // LCP length.
         index[l] + lcpLength >= targetLength or  // This shouldn't
                                                     // happen
         // End on a mismatch.
         ThreeBit[query[lcpLength]] >= 4 or 
				 Compare::Compare(target[index[l] + lcpLength], query[lcpLength]) != 0
         
         ) {
			 break;
		 }
     

     
		 //
		 // Store the bounds for the current offset.  These are used later
		 // to expand the search if necessary.
		 //
		 lcpLeftBounds.push_back(l);
		 lcpRightBounds.push_back(r);
     lcpLength++;


	 }
	 return lcpLength;
 }


 int SearchLow(T *target, T *query, DNALength queryLength, SAIndex l, SAIndex r, SAIndex &low, unsigned int offset=0) {

	 long midPos;
	 int high;
	 int numSteps = 0;
	 // 
	 // Boundary conditions, the string is either before (lexicographically) the text
	 // or after.
	 //
	 if (StringLessThanEqual(&query[offset], 
													 queryLength-offset, 
													 &target[index[l]+offset], 
													 length - index[l]-offset)) {
		 low = l;
		 return low;
	 }
	 else if (StringLessThan(&target[index[r]+offset], 
													 length - index[r] - offset , 
													 &query[offset], 
													 queryLength  - offset)) {
		 low = length;
		 return low;
	 }

	 //
	 // The string fits somewhere in the text.
	 //
	 low = l;
	 high = r;
	 long diff = ((long) high) - ((long) low);
	 while (diff > 1) {
		 ++numSteps;
		 midPos = ((long) high) + ((long) low);
		 midPos = midPos / 2;
		 if (StringLessThanEqual(&query[offset], queryLength-offset, &target[index[midPos]+offset], length - index[midPos]-offset)) {
			 high = midPos;
		 }
		 else {
			 low = midPos;
		 }
		 diff = ((long) high) - ((long) low);
	 }

	 //
	 // The search is for the least position such that the query is greater than or equal to the text.
	 // High tracks the positions that may be equal to the query, and low is strictly less than the query. 
	 // At the end of the search, high is either pointing to the query, or the first element where the query
	 // could be placed before high without changing the order of target.
	 //
	 low = high;
	 diff = ((long) high) - ((long) low);
	 return low;
	 //		cout << "search low took: " << numSteps << endl;
 }


 int SearchHigh(T *target, T *query, DNALength queryLength, SAIndex l, SAIndex r,  SAIndex &high, unsigned int offset=0) {

	 //
	 // Find the last position where the query is less than the target.
	 //
	 long midPos;
	 int low;
	 int numSteps = 0;
	 // 
	 // Boundary conditions, the string is either before (lexicographically) the text
	 // or after.
	 //
	 if (StringLessThan(&target[index[r]+offset], length - index[r] - offset, &query[offset], queryLength-offset)) {
		 high = -1;
		 return high;
	 }

	 //
	 // The string fits somewhere in the text.
	 //
	 low = l;
	 high = r;
	 long diff = ((long) high) - ((long) low);
	 while (diff > 1) {
		 ++numSteps;
		 midPos = ((long) high) + ((long) low);
		 midPos = midPos / 2;
		 if (StringLessThan(&query[offset], queryLength - offset, &target[index[midPos]+offset], length - index[midPos] - offset)) {
			 high = midPos;
		 }
		 else {
			 low = midPos;
		 }
		 diff = ((long) high) - ((long) low);
	 }

	 //
	 // The search is for the last position where the query is less than or equal to the text.  High is 
	 // strictly greater than or the query.  Low is less than or equal to the query.  At the end, low will be 
	 // the last spot where query could be inserted after and not wreck the ordering of the array.
	 //
	 high = low;
	 //		cout << "search high took: " << numSteps << " steps." << endl;
 }
};



#endif
