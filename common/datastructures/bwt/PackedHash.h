#ifndef DATASTRUCTURES_BWT_PACKED_HASH_H_
#define DATASTRUCTURES_BWT_PACKED_HASH_H_

#include <vector>
#include <map>
#include <iostream>
#include "../../Types.h"
#include "../../utils.h"
#include "../../utils/BitUtils.h"
#include "../../DNASequence.h"

class PackedHash {
 public:

	int tableLength;
	uint32_t *table;
	uint64_t *values;
	vector<int> hashLengths;
	static const uint32_t BinNumBits = 5;
	static const uint32_t BinSize = 1 <<(BinNumBits);
	
	/*
	 * Create a mask that retains the lower 5 bits (0 .. 31) of a
	 * position so that pos % 32 may be computed by a shift.
	 */
	static const uint32_t BinModMask = 0x1FU; 
	
	void Allocate(uint32_t sequenceLength) {
		tableLength = CeilOfFraction(sequenceLength, (DNALength) BinSize);
		table  = new uint32_t[tableLength];
		values = new uint64_t[tableLength];
		std::fill(&table[0], &table[tableLength], 0);
		std::fill(&values[0], &values[tableLength], 0);
		hashLengths.resize(tableLength);
		std::fill(hashLengths.begin(), hashLengths.end(), 0);
	}

	void PrintBinSummary() {
		//
		// Report some stats on the hash.
		//
		DNALength p;
		map<int,int> countMap;
		int card;
		for (p = 0; p < tableLength; p++ ){ 
			card = CountBits(table[p]);
			countMap[card]++;
		}
		map<int,int>::iterator mapit;
		for (mapit = countMap.begin(); mapit != countMap.end(); ++mapit) {
			cout << mapit->first << " " << mapit->second << endl;
		}
	}

	uint32_t LookupBinAtPos(DNALength pos) {
		/* 
		 * Each bucket contains 32 positions.  Membership is simply when
		 * the bit at pos & BinModMask is set.  There should never be
		 * collisions of multiple positions (from different areas in the
		 * genome) mapping to the same position in a bin.
		 */
		return table[pos/BinSize] & (1 << (pos & BinModMask));
	}
	
	void ValueToList(uint64_t &storage, DNALength newValue, int newValuePos) {

		/*
		 * This is called when one is attempting to add a spot to storage,
		 * but there are already two values stored in it, so storage must
		 * be converted to a pointer, and the values added to a list on
		 * the heap.  The size of the list is necessarily 3 at this point,
		 * because there are two values in storage that must be moved to
		 * the list, and the new value as well. 
		 * The values are copied to the list in sorted order, and sorting
		 * is handled case-by-case since there are only 3 cases.
		 */

		DNALength v0, v1;
		v0 = ((DNALength)storage);
		v1 = ((DNALength)(storage >> 32));
		DNALength *storagePtr = new DNALength[3];
		storage = (uint64_t) storagePtr;
		
		//
		// Only a couple of options, so handle them directly here
		//
		if (newValuePos == 0)	{
			storagePtr[0] = newValue; 
			storagePtr[1] = v0; storagePtr[2] = v1;
		}
		else if (newValuePos == 1) {
			storagePtr[0] = v0; storagePtr[1] = newValue; storagePtr[2] = v1;
		}
		else if (newValuePos == 2) {
			storagePtr[0] = v0; storagePtr[1] = v1; storagePtr[2] = newValue;
		}
		else {
			assert("ERROR! Somehow expected to only add 3 elements to an array of length 3, but the position of the new value is greater than the length of the array" && 0);
		}
	}

	void InsertValueInList(uint64_t &storage, int curStorageLength, DNALength newValue, int newValuePos) {
		/*
		 * This simply creates a new list with size 1 larger than before,
		 * and inserts the new value into its position that maintains
		 * sorted order in the list.
		 */
		DNALength *newListPtr = new DNALength[curStorageLength + 1];
		//
		// Copy the values from the old list making space for the new
		// value.
		// 
		if (newValuePos > 0) {
			memcpy(newListPtr, ((DNALength*)storage), newValuePos * sizeof(DNALength));
		}
		if (newValuePos < curStorageLength) {
			memcpy(&newListPtr[newValuePos+1], &((DNALength*)storage)[newValuePos], (curStorageLength - newValuePos)*sizeof(uint32_t));
		}
		assert(curStorageLength < 32);
		newListPtr[newValuePos] = newValue;
		delete[] ((DNALength*)storage);
		storage = (uint64_t)newListPtr;
	}
		
	void DirectlyStoreValue(uint64_t &storage, int curStorageLength, DNALength value, int valuePos) {
		/*
		 * In this instance, the value may be copied to either the first
		 * half or the second half of 'storage', without having to turn
		 * storage into a list.  The values must be stored in sorted
		 * order.  If 'storage' is empty, the correct place is at the
		 * beginning of storage.  If there already is a value at the
		 * beginning, it may be necessary to shift the existing value over
		 * to keep everything in sorted order. 
		 */
		
		if (curStorageLength == 0) {
			//
			// Nothing here, just store the value.
			//
			storage = value;
		}
		else if (valuePos == 0) {
			//
			// Place the value at the beginning of the storage.
			//
			storage = storage << 32;
			storage = storage + value;
		}
		else {
			// 
			// Place the value at the end of storage.
			//
			uint64_t longValue = value;
			longValue = longValue << 32;
			storage += longValue;
		}
	}

	int AddValue(DNALength pos, DNALength value) {
		//
		// The bucket is either a values[pos] that can store up to two
		// values, or it is a pointer to a list of values.  The values are
		// always in order of their corresponding position in the BWT
		// string. 
		// To add a value to this bucket, first check to see if
		// values[pos] has enough room to simply put it there, otherwise,
		// either a list already exists and the value must be inserted, or
		// the values[pos] must be converted to a list, and then the new
		// value inserted.

		//
		// First, operate on the assumption that no value gets added
		// twice.
		//
		UInt bin = pos / BinSize;
		UInt bit = pos & BinModMask;
		assert((table[bin] & ( 1 << bit)) == 0);
		
		//
		// Now, add the pos and determine where it is in the array.
		//
		table[bin] = table[bin] + ( 1 << bit);

		//
		// Mask off everything above this bit.
		//
		UInt mask = (UInt)-1;
		mask >>= (31 - bit);
		UInt lowerBits = table[bin] & mask;
		int  rank = CountBits(lowerBits) - 1;
		int  card = CountBits(table[bin]);

		if (card < 3) {
			DirectlyStoreValue(values[bin], card-1, value, rank);
		}
		else if (card == 3) {
			ValueToList(values[bin], value, rank);
		}
		else {
			InsertValueInList(values[bin], card-1, value, rank);
		}
		return card;
	}

	int LookupValue(DNALength pos, DNALength &value) {

		/* 
		 * Check to see if there is a value stored for 'pos'.  If so,
		 * store it in value, and return 1 for success.
		 */
		UInt binIndex = pos/BinSize;
		UInt bin      = table[binIndex];
		UInt setBit = bin & (1 << (pos & BinModMask));
		if (setBit == 0) {
			return 0;
		}

		// 
		// GetSetBitPosition64 returns the position relative to the most
		// significant bit in a 64-bit word, starting at MSB = 1.  This
		// should never be less than 32, since it's counting bits in a 32
		// bit word.  
		//

		int  bitPos     = GetSetBitPosition32(setBit);
		UInt bitPosMask = ((UInt)-1) >> (32 - bitPos-1);;
		
		//
		// Determine how to interpret the value of this bucket.  It is
		// either a pointer or two words.  If there are more than two bits
		// set in this bucket, it is a pointer.  Otherwise, pick the half
		// of the 64 bit word based on the index in the bucket.
		//
		int nSet         = CountBits(bin);
		int bitRank      = CountBits(bin & bitPosMask) - 1;
		assert(nSet > 0);
		if (nSet <= 2) {
			// return lower 32 bits
			if (bitRank == 0) {
				value = ((uint32_t)values[binIndex]);
				return 1;
			}
			else {
				value = ((uint32_t)(values[binIndex]>>32));
				return 1;
			}
		}
		else {
			/*
			 * In this instance, values is a pointer rather than a pair of values.
			 */
			value = (uint32_t) ((DNALength*)values[binIndex])[bitRank];
			return 1;
		}
		//
		// This is reached if nothing is stored in value. For now, this
		// shouldn't be reached, so make sure the program bails here.
		//
		assert(0);
		return 0;
	}	

	/*
	 * Define binary I/O routines for the packed hash. The sizes of the
	 * tables are constant, but then the values table has some lists.
	 * Those are written past the end of the tables treating the last
	 * half of the file as a heap structure.  They are simlarly read in.
	 *
	 */

	void Write(ostream &out) {
		out.write((char*) &tableLength,sizeof(tableLength));
		if (tableLength > 0) {
			out.write((char*) table, tableLength * sizeof(table[0]));
			out.write((char*) values, tableLength * sizeof(values[0]));
		}
		int tablePos;
		for (tablePos = 0; tablePos < tableLength; tablePos++) {
			int nSetBits = CountBits(table[tablePos]);
			if( nSetBits > 2) {
				out.write((char*)values[tablePos], sizeof(uint32_t)*nSetBits);
			}
		}
	}

	void Read(istream &in) {
		in.read((char*)&tableLength, sizeof(tableLength));
		if (tableLength > 0) {
			table  = new uint32_t[tableLength];
			values = new uint64_t[tableLength];
			in.read((char*)table, sizeof(uint32_t)*tableLength);
			in.read((char*)values, sizeof(uint64_t)*tableLength);
			int tablePos;
			for (tablePos = 0; tablePos < tableLength; tablePos++) {
				int nSetBits = CountBits(table[tablePos]);
				if (nSetBits > 2) {
					values[tablePos] = (uint64_t) new uint32_t[nSetBits];
					in.read((char*)values[tablePos], nSetBits * sizeof(uint32_t));
				}
			}
		}
	}
};

#endif
