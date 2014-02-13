#ifndef DATASTRUCTURES_READS_READ_LIST_H_
#define DATASTRUCTURES_READS_READ_LIST_H_
#include <vector>
#include <string>
#include "../../Types.h"

template<typename T_Read>
class CompareReadTitle {
 public:
	int operator()( T_Read lhs, const string &rhs) const {
		return lhs->GetTitle().compare(rhs) < 0;
  }
};


template<typename T_Read>
class CompareReadsByTitles {
public:
  int operator()(const T_Read *lhs, const T_Read* rhs) const {
	  return lhs->GetTitle().compare(rhs->GetTitle()) < 0;
  }
};
			
template<typename T_Read>
class ReadList {
public:
	vector<T_Read*> reads;
	typedef vector<T_Read*> ReadPtrList;
	int LookupReadByName(string readName, VectorIndex &index) {
		typename ReadPtrList::iterator readIndexIt;
		readIndexIt = std::lower_bound(reads.begin(), reads.end(), readName, CompareReadTitle<T_Read*>());
		if (readIndexIt == reads.end()) return 0;
		else {
			index = readIndexIt - reads.begin();
			if ((*readIndexIt)->GetTitle().compare(readName) != 0) { return 0; } 
			else { return 1; }
		}
	}
	void Add(const T_Read *readPtr) {
		reads.push_back((T_Read*) readPtr);
	}
	void Order() {
		std::sort(reads.begin(), reads.end(), CompareReadsByTitles<T_Read>());
	}
};

#endif
