#ifndef HASHED_TUPLE_LIST_H_
#define HASHED_TUPLE_LIST_H_
#include "TupleList.h"
#include "TupleMetrics.h"
#include "../DNASequence.h"

template<typename T_Tuple>
class HashedTupleList {
 public:
	long mask;
	vector<TupleList<T_Tuple> > hashTable;
	int hashLength;
  int hashTableLength;
  typedef T_Tuple Tuple;
  //
  // Provide a defalt constructor with a small tuple size for testing.
  //
  HashedTupleList() {
    Initialize(5);
  }

  void Initialize(int _hashLength) {
		mask = 0;
		int i;
		hashLength = _hashLength;
		hashTable.resize(1L << (hashLength*2));
    hashTableLength = hashTable.size();

		for (i = 0; i < hashLength; i++ ){
			mask = mask << 2L;
			mask = mask + 3L;
		}
  }
      
	HashedTupleList(int _hashLength) {
    Initialize(_hashLength);
    cout << hashTable.size() << endl;
	}

  void clear() {
    // Synonym.
    Clear();
  }

	void Clear() {
		int i;
		for (i = 0; i < hashTableLength; i++ ){ 
			hashTable[i].tupleList.clear();
		}
	}

  void Sort() {
		int i;
		for (i = 0; i < hashTableLength; i++ ){ 
			sort(hashTable[i].tupleList.begin(), hashTable[i].tupleList.end());
		}
  }    

  void Append(T_Tuple tuple) {
		int hashValue = tuple.tuple & mask;
    cout << "htl adding " << tuple.tuple << endl;
		hashTable[hashValue].tupleList.push_back(tuple);
	}
    
	void Insert(T_Tuple tuple) {
    cout << "htl adding " << tuple.tuple << endl;
		int hashValue = tuple.tuple & mask;
		hashTable[hashValue].Insert(tuple);
	}

  int Find(T_Tuple tuple) {
    int hashValue, index;
    return Find(tuple, hashValue, index);
  }

  void Print() {
    int i;
		for (i = 0; i < hashTableLength; i++ ){ 
      hashTable[i].Print();
    }
  }
  //
  // Provide a version of find that stores easy access to the original
  // tuple.
  //
	int Find(T_Tuple tuple, int &hashValue, int &index) {
		hashValue = tuple.tuple & mask;
		if (hashTable[hashValue].size()) {
      return ((index = hashTable[hashValue].Find(tuple)) != -1);
		}
		else {
			return 0;
		}
	}

  void FindAll(T_Tuple &tuple, 
               typename vector<T_Tuple>::const_iterator &firstPos, 
               typename vector<T_Tuple>::const_iterator &endPos ) {
    int hashValue;
    hashValue = tuple.tuple & mask;
    hashTable[hashValue].FindAll(tuple, firstPos, endPos);
  }


  int GetHashLength() {
    return hashLength;
  }
};


template<typename T_Tuple>
void SequenceToHash(DNASequence &seq, HashedTupleList<T_Tuple> &hash, TupleMetrics &tm) {
	int i;
	T_Tuple tuple;
  int res = 0;
	for (i = 0; i < seq.length - hash.hashLength + 1; i++ ) {
    if ((res and (res = tuple.ShiftAddRL(seq.seq[i+tm.tupleSize-1], tm))) or
        (!res and (res = tuple.FromStringRL(&seq.seq[i], tm)))) {
      hash.Insert(tuple);
    }
  }
}




#endif
