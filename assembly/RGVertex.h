#ifndef ASSEMBLY_RG_VERTEX_H_
#define ASSEMBLY_RG_VERTEX_H_

#include <vector>
#include <iostream>
#include "Types.h"
using namespace std;


template<typename T_Key>
class RGVertex {
	public:
	T_Key key;
	vector<int> outEdges;
	vector<int> inEdges;
	bool traversed;
	int flowIn;
	int flowOut;
	int flowBalance;
	int isRepeat;
	RGVertex<T_Key>() {
		traversed = false;
		flowIn = flowOut = 0;
		isRepeat = false;
	}
	RGVertex<T_Key>(T_Key &keyP) {
		key = keyP;
		traversed = false;
		isRepeat = false;
		flowIn = flowOut = 0;
	}
	int InDegree() {
		return inEdges.size();
	}
	int OutDegree() {
		return outEdges.size();
	}
  int IsCross() {
		return InDegree() == 2 and OutDegree() == 2;
  }

	RGVertex<T_Key>&operator=(const RGVertex<T_Key> &rhs) {
		key = rhs.key;
		traversed = rhs.traversed;
		inEdges = rhs.inEdges;
		outEdges = rhs.outEdges;
		return *this;
	}
	RGVertex<T_Key>(const RGVertex<T_Key> &rhs) {
		*this = rhs;
	}
			
	friend ostream& operator<<(ostream &out, const RGVertex<T_Key> &rhs) {
		out << rhs.key << " ";
		VectorIndex inIndex;
		out << rhs.inEdges.size() << " ";
		for (inIndex = 0; inIndex < rhs.inEdges.size(); inIndex++) {
			out << rhs.inEdges[inIndex] << " ";
		}

		VectorIndex outIndex;
		out << rhs.outEdges.size() << " ";
		for (outIndex = 0; outIndex < rhs.outEdges.size(); outIndex++) {
			out << rhs.outEdges[outIndex] << " ";
		}
		out << endl;
		return out;
	}
 
	friend istream & operator>>(istream &in, RGVertex< T_Key > &rhs) {
		
		in >> rhs.key;
		VectorIndex inIndex;
		VectorIndex numInEdges;
		in >> numInEdges;
		if (numInEdges > 0) {
			rhs.inEdges.resize(numInEdges);
			for (inIndex = 0; inIndex < rhs.inEdges.size(); inIndex++) {
				in >> rhs.inEdges[inIndex];
			}
		}
		VectorIndex outIndex;
		VectorIndex numOutEdges;
		in >>numOutEdges;
		if (numOutEdges > 0) {
			rhs.outEdges.resize(numOutEdges);
			for (outIndex = 0; outIndex < numOutEdges; outIndex++) {
				in >> rhs.outEdges[outIndex];
			}
		}
		return in;
 }		
	void SetKey(T_Key &keyP) {
		key = keyP;
	}
};
	

#endif
