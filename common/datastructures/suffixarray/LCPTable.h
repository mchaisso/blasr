#ifndef LCP_TABLE_H_
#define LCP_TABLE_H_


#include <map>

using namespace std;
template <typename T>
class LCPTable {
	//
	// Change the following TWO type defs if 
	// the max LCP is changed.
	//
	typedef short SignedPrefixLength;
	typedef unsigned short PrefixLength;
	typedef map<int, int> LongPrefixMap;
	PrefixLength maxPrefixLength;
	LongPrefixMap llongPrefixMap, rlongPrefixMap;
	int tableLength;
 public:
	PrefixLength *llcp, *rlcp;
	LCPTable() {
		tableLength =0;
		llcp = rlcp = NULL;
	}

	LCPTable(T* data, unsigned int pTableLength) {
		Init(data, pTableLength);
	}

	inline	int LengthLongestCommonPrefix(T* a, int alen, T* b, int blen) {
		int i;
		for (i = 0 ; i < alen and i < blen; i++ ) 
			if (a[i] != b[i])
				break;
		return i;
	}

	void Init(T* data, unsigned int pTableLength, unsigned int *index) {
		tableLength = pTableLength;
		maxPrefixLength = (PrefixLength) (SignedPrefixLength(-1));
		llcp = new PrefixLength[tableLength];
		rlcp = new PrefixLength[tableLength];
		std::fill(llcp, llcp + tableLength, 0);
		std::fill(rlcp, rlcp + tableLength, 0);
		FillTable(data, index);
	}
	
	int SetL(int index, int length) {
		assert(index >= 0);
		assert(index < tableLength);
		if (index >= maxPrefixLength) {
			llcp[index] = maxPrefixLength;
			llongPrefixMap[index] = length;
		}
		else {
			llcp[index] = length;
		}
		return llcp[index];
	}

	int SetR(int index, int length) {
		assert(index >= 0);
		assert(index < tableLength);
		if (index >= maxPrefixLength) {
			rlcp[index] = maxPrefixLength;
			rlongPrefixMap[index] = length;
		}
		else {
			rlcp[index] = length;
		}
		return rlcp[index];
	}


	int GetL(int index) {
		if (llcp[index] == maxPrefixLength) {
			assert(llongPrefixMap.find(index) != llongPrefixMap.end());
			return llongPrefixMap[index];
		}
		else {
			return llcp[index];
		}
	}

	int GetR(int index) {
		if (rlcp[index] == maxPrefixLength) {
			assert(rlongPrefixMap.find(index) != llongPrefixMap.end());
			return rlongPrefixMap[index];
		}
		else {
			return rlcp[index];
		}
	}


	~LCPTable() {
		/*
			if (llcp != NULL) 
			delete[] llcp;
			llcp = NULL;
			if (rlcp != NULL)
			delete[] rlcp;
			rlcp = NULL;
		*/
		// the two maps automatically go away.
	}

	void WriteLCPTable(ofstream &out) {
		out.write((char*) &tableLength, sizeof(tableLength));
		out.write((char*) llcp, sizeof(PrefixLength)*tableLength);
		out.write((char*) rlcp, sizeof(PrefixLength)*tableLength);
		typename LongPrefixMap::iterator llcpTableIt, llcpTableEnd;
		int llongPrefixMapSize = llongPrefixMap.size();
		out.write((char*) &llongPrefixMapSize, sizeof(llongPrefixMapSize));
		for((llcpTableIt = llongPrefixMap.begin(),
				 llcpTableEnd = llongPrefixMap.end());
				llcpTableIt != llcpTableEnd;
				++llcpTableIt) {
			out.write((char*) &llcpTableIt->first, sizeof(typename LongPrefixMap::key_type));
			out.write((char*) &llcpTableIt->second, sizeof(typename LongPrefixMap::mapped_type));
		}
		typename LongPrefixMap::iterator rlcpTableIt, rlcpTableEnd;
		int rlongPrefixMapSize = rlongPrefixMap.size();
		out.write((char*) &rlongPrefixMapSize, sizeof(rlongPrefixMapSize));
		for((rlcpTableIt = rlongPrefixMap.begin(),
				 rlcpTableEnd = rlongPrefixMap.end());
				rlcpTableIt != rlcpTableEnd;
				++rlcpTableIt) {
			out.write((char*) &rlcpTableIt->first, sizeof(typename LongPrefixMap::key_type));
			out.write((char*) &rlcpTableIt->second, sizeof(typename LongPrefixMap::mapped_type));
		}
	}

	void ReadLCPTables(ifstream &in) {
		in.read((char*) &tableLength, sizeof(tableLength));
		in.read((char*) &llcp, sizeof(PrefixLength)*tableLength);
		in.read((char*) &rlcp, sizeof(PrefixLength)*tableLength);
		int longPrefixMapSize;
		in.read((char*) &longPrefixMapSize, sizeof(longPrefixMapSize));
		int i;
		int index, lcpLength;
		// The rest is stored as a tree, but read and construct this on the fly.
		for (i = 0; i < longPrefixMapSize; i++) {
			in.read((char*) &index, sizeof(index));
			in.read((char*) &lcpLength, sizeof(lcpLength));
			llongPrefixMap[index] = lcpLength;
		}

		// Read the rlcp
		in.read((char*) &longPrefixMapSize, sizeof(longPrefixMapSize));
		for (i = 0; i < longPrefixMapSize; i++) {
			in.read((char*) &index, sizeof(index));
			in.read((char*) &lcpLength, sizeof(lcpLength));
			rlongPrefixMap[index] = lcpLength;
		}
	}


	void FillTable(T* data, unsigned int *index) {
		//
		// This assumes that the index table is now in sorted order.
		//
		FillTable(0, tableLength, data, index);
	}

	void FillTable(unsigned int low, unsigned int high, T* data, unsigned int *index) {
		if (low == high) 
			return;
		unsigned int mid = (low + high) / 2;
		assert (mid != low);
		SetL(mid, LengthLongestCommonPrefix(&data[index[low]], tableLength - index[low],
																				&data[index[mid]], tableLength - index[mid]));
		FillTable(low, mid, data, index);
		assert(mid != high);
		SetR(mid, LengthLongestCommonPrefix(&data[index[mid]], tableLength - index[mid],
																				&data[index[high]], tableLength - index[high]));
		FillTable(mid+1, high, data, index);
	}

};

#endif
