#ifndef TUPLE_LIST_H_
#define TUPLE_LIST_H_

#include <string>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "../Types.h"
#include "TupleMetrics.h"

using namespace std;

template<typename T>
class TupleList {
  int listLength;
	TupleMetrics tm;
 public:
  typedef T Tuple;
  vector<T> tupleList;

  void Reset() {
    vector<T>().swap(tupleList);
  }

  TupleList() {
		listLength = 0;
  }
  T &operator[](int index) {
		return tupleList[index];
  }

	void GetTupleMetrics(TupleMetrics &ptm) {
		ptm = tm;
	}
	void SetTupleMetrics(TupleMetrics &ptm) {
		tm = ptm;
	}
	int size() {
		return tupleList.size();
	}
  int GetLength() {
		return tupleList.size();
  }
  int InitFromFile(string &fileName) {
		ifstream listIn;
		listIn.open(fileName.c_str(), ios_base::binary);
		if (!listIn)
			return 0;
		listIn.read((char*) &listLength, sizeof(int));
		listIn.read((char*) &tm.tupleSize, sizeof(int));
		tm.InitializeMask();
		//list = new T[listLength];
		tupleList.resize(listLength);
		listIn.read((char*) &tupleList[0], sizeof(T) * listLength);
		return 1;
  }

	void clear() {
		tupleList.clear();
		listLength = 0;
	}
	
  int WriteToFile(string &fileName) {
		ofstream listOut;
		listOut.open(fileName.c_str(), ios_base::binary);
		if (!listOut)
			return 0;
		listLength = tupleList.size();
		cout << "writing tuple lis of length " << listLength << endl;
		listOut.write((char*) &listLength, sizeof(int));
		listOut.write((char*) &tm.tupleSize, sizeof(int));
		listOut.write((char*) &tupleList[0], sizeof(T)*listLength);
		return 1;
  }

  //
  // Find one instance of a match.
  //
  int Find( T& tuple)  {
    typename vector<T>::const_iterator begin, end, matchIt;
    begin = tupleList.begin();
    end   = tupleList.end();
    matchIt = lower_bound(begin, end, tuple);
    if (*matchIt != tuple) {
      return -1;
    }
    else {
      return matchIt - tupleList.begin();
    }
	}

  //
  // Find the boundaries of all instances of a match.
  //
  void FindAll(T &tuple, typename vector<T>::const_iterator &firstPos, typename vector<T>::const_iterator &endPos ) {
    firstPos = lower_bound(tupleList.begin(), tupleList.end(), tuple);
    typename vector<T>::const_iterator firstPos2;
    endPos = tupleList.end();
    endPos = upper_bound(firstPos, endPos, tuple);
    while (endPos != tupleList.end()) {
      if (*endPos != tuple) {
        return;
      }
      else {
        endPos++;
      }
    }
  }

	void Append( T&tuple) {
		tupleList.push_back(tuple);
	}
	
	void Insert(T&tuple) {
		// insert and maintain order.
		typename vector<T>::iterator pos;
		pos = std::lower_bound(tupleList.begin(), tupleList.end(), tuple);
		tupleList.insert(pos, tuple);
	}

	void Sort() {
		sort(tupleList.begin(), tupleList.end());
	}
	

  void Print() {
    int i;
    for (i = 0; i< tupleList.size(); i++) {
      cout << tupleList[i].tuple << endl;
    }
  }
};




#endif
