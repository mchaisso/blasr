#ifndef ASSEMBLY_PATH_COLLECTION_H_
#define ASSEMBLY_PATH_COLLECTION_H_

#include "Path.h"
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>


using namespace std;
template<typename T_Key>
class PathCollection {
public:
 vector<Path<T_Key> > paths;
 void AddPath(Path<T_Key> &path) {
		paths.push_back(path);
 }
   int size() { return paths.size();} 
 void Condense() {
	 sort(paths.begin(), paths.end(), PathLessThan<T_Key>());
	 UInt curIndex, packedIndex;
   curIndex = 0;
   packedIndex = 0;
   while (curIndex < paths.size()) {
		 UInt equalPathIndex = curIndex;
     while (equalPathIndex < paths.size() and paths[equalPathIndex] == paths[curIndex]) {	
		    equalPathIndex++;
     }
     paths[curIndex].multiplicity = equalPathIndex - curIndex;
		 paths[packedIndex] = paths[curIndex];
		 curIndex = equalPathIndex;
	   packedIndex++;
   }
   paths.resize(packedIndex);	 
 }
 void RemoveShortPaths(int minLength) {
		UInt curIndex, packedIndex;	
    for(curIndex=0,packedIndex = 0;
				curIndex < paths.size(); curIndex++ ) {
			if (paths[curIndex].size() >= minLength) {
				paths[packedIndex] = paths[curIndex];
				packedIndex++;
			}
		}
		paths.resize(packedIndex);
  }

  void Read(istream &in, int minCount = 2) {
    string line;
    while(getline(in,line)) {
      stringstream linestrm;
			linestrm.str(line);
		  Path<T_Key> path;
			linestrm >> path.multiplicity;
			if (path.multiplicity < minCount) { continue; }
			int e;
      while((linestrm >> e)) {
				path.edges.push_back(e);
			}
      paths.push_back(path);
     }
  }

 void Print(ostream &out) {
		UInt p;
    for (p = 0; p < paths.size(); p++ ) {
			out << paths[p].multiplicity << " ";
      UInt e;
      for (e = 0; e < paths.size(); e++ ){ 
				out << paths[p].edges[e] << " ";
      }
      out << endl;
    }
 }
 

 void RemoveEnclosedPaths() {
   vector<bool> enclosed;
	 enclosed.resize(paths.size());
	 fill(enclosed.begin(), enclosed.end(), false);
   Path<T_Key> subpath;
   UInt p;
   for (p = 0; p < paths.size(); p++ ) {
     // count all the paths contined by this one.
     UInt pe;
     for (pe = 0; pe < paths[p].size();   pe++) {
			 

     }
   } 

 }

};



#endif
