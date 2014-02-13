#ifndef ASSEMBLY_PATH_H_
#define ASSEMBLY_PATH_H_
#include <vector>
#include "Types.h"
using namespace std;
template<typename T_Key>
class Path {
public:
	vector<T_Key> edges;
	int multiplicity;
  int size() const { return edges.size(); }
	UInt operator[](UInt index) {
		return edges[index];
  }
  Path& operator=(const Path &rhs) {
	  edges = rhs.edges;
	  multiplicity = rhs.multiplicity;
    return *this;
  }

  int operator==(const Path &rhs) const {
		UInt pos;
		if (rhs.size() != edges.size()) { return 0; }
		for (pos = 0; pos < rhs.size() && pos < edges.size(); pos++) {
	    if (edges[pos] != rhs.edges[pos]) return 0;
    }
		return 1;
  }
	void AddEdge(UInt edgeIndex) {
		edges.push_back(edgeIndex);
  }

};
template<typename T_Key>
class PathLessThan { 
public:
	int operator()(const Path<T_Key> &lhs, const Path<T_Key> &rhs) const {
		UInt pos;
		for (pos = 0; pos < lhs.size() and pos < rhs.size(); pos++) {
			if (lhs.edges[pos] < rhs.edges[pos]) { return 1; }
		  if (lhs.edges[pos] > rhs.edges[pos]) { return 0; } 
    }
		if (lhs.size() < rhs.size() ) { return 1;}
		// if here, rhs.size() <= lhs.size(), and the.edgess are otherwise equal, 
		// so lhs > rhs.size
    return 0;
  }
};

#endif
