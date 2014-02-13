#ifndef SHARED_SUFFIX_ARRAY_H_
#define SHARED_SUFFIX_ARRAY_H_

#include <stdlib.h>
#include <sstream>
#include "SuffixArray.h"
#include "../../ipc/SharedMemoryAllocator.h"
#include "../../tuples/DNATuple.h"
#include "../../tuples/CompressedDNATuple.h"
#include "../../algorithms/compare/CompareStrings.h"


template<typename T, 
	typename Sigma,
	typename Compare = DefaultCompareStrings<T>,
	typename Tuple   = DNATuple >
class SharedSuffixArray : public SuffixArray<T, Sigma, Compare, Tuple> {
	string shmIdTag;
 public:
	SAIndex *indexShared;
	int indexID;
	string indexHandle;
	SAIndex *lookupTableShared;
	int lookupTableID;
	string lookupTableHandle;
	int lookupPrefixLength;
	
	void InitShmIdTag() {
		stringstream tagStrm;
		tagStrm << "_" << getpid();
		shmIdTag = tagStrm.str();
	}
	

	void ReadSharedArray(ifstream &saIn) {
		cout << "reading a shared suffix array index." << endl;
	 saIn.read((char*) &this->length, sizeof(int));
	 indexHandle = "suffixarray.index." + shmIdTag;
	 AllocateMappedShare(indexHandle, this->length + 1, indexShared, indexID);
		cout << "the shared index is: " << indexShared << endl;
	 this->index = indexShared;
		cout << "the index used is: " << this->index << endl;
	 this->ReadAllocatedArray(saIn);
	}

	void ReadSharedLookupTable(ifstream &saIn) {
		this->ReadLookupTableLengths(saIn);
		lookupTableHandle = "suffixarray.lookuptable." + shmIdTag;
		AllocateMappedShare(lookupTableHandle, this->lookupTableLength + 1, lookupTableShared, lookupTableID);
		this->lookupTable = lookupTableShared;
		this->ReadAllocatedLookupTable(saIn);
	}

	void ReadShared(string &inFileName) {
	 ifstream saIn;
	 InitShmIdTag();
	 saIn.open(inFileName.c_str(), ios::binary);
	 this->ReadComponentList(saIn);
	 if (this->componentList[SuffixArray<T,Sigma,Compare,Tuple>::CompArray]) {
		 this->ReadSharedArray(saIn);
	 }
	 if (this->componentList[SuffixArray<T,Sigma,Compare,Tuple>::CompLookupTable]) {
		 this->ReadSharedLookupTable(saIn);
	 }
	 saIn.close();
 }

	void FreeShared() {

	 if (this->componentList[SuffixArray<T,Sigma,Compare,Tuple>::CompArray]) {
		 shm_unlink(indexHandle.c_str());
	 }
	 if (this->componentList[SuffixArray<T,Sigma,Compare,Tuple>::CompLookupTable]) {
		 shm_unlink(lookupTableHandle.c_str());
	 }
	}


};


#endif
