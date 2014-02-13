#ifndef ASSEMBLY_RG_EDGE
#define ASSEMBLY_RG_EDGE

using namespace std;
#include "Types.h"

class RGEdge {
 public:
	VectorIndex src;
	VectorIndex dest;
	bool connected;
  bool traversed;
	int count;
	int flow;
	int minFlow;
	RGEdge() {
		connected = false;
		count = 0;
		flow  = 0;
		minFlow = 0;

	}

	RGEdge(const VectorIndex srcP, const VectorIndex destP){ 
		src  = srcP;
		dest = destP;
		count = 0;
		flow  = 0;
		minFlow = 0;
	}

	RGEdge&operator=(const RGEdge &rhs) {
		src = rhs.src;
		dest = rhs.dest;
		count = rhs.count;
		flow = rhs.flow;
		minFlow = rhs.minFlow;
		return *this;
	}

	RGEdge(const RGEdge &rhs) {
		(*this) = rhs;
	}
	
	friend ostream &operator<<(ostream &out, const RGEdge &rhs) {
		out << rhs.src << " " << rhs.dest << " " << rhs.count << endl;
		return out;
	}
	
	friend istream &operator>>(istream &in, RGEdge &rhs) {
		in >> rhs.src >> rhs.dest >> rhs.count;
		return in;
	}
};



#endif
