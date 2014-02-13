#include "RepeatGraph.h"
#include "RGEdgeOperations.h"
#include <string>
#include <set>
#include <iostream>
#include <sstream>
using namespace std;

void RecursivePrintComponent(RepeatGraph<string> &rg, int curVertex, ostream &out) {
  int outEdgeIndex, outEdge, inEdgeIndex, inEdge;
  int dest, src;
  rg.vertices[curVertex].traversed = true;
cout << curVertex << endl;
  for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[curVertex].outEdges.size(); outEdgeIndex++) {
		outEdge =rg.vertices[curVertex].outEdges[outEdgeIndex];
    dest = rg.edges[outEdge].dest;
    if (rg.edges[outEdge].traversed == false) {
			out << "\"" << rg.vertices[curVertex].key << "\" -> \"" << rg.vertices[dest].key << "\" [label=\""<< rg.edges[outEdge].count << "\"];" << endl;
			rg.edges[outEdge].traversed = true;
      if (rg.vertices[dest].traversed == false) {
        RecursivePrintComponent(rg, dest, out);
      }
    }
  }
  for (inEdgeIndex = 0; inEdgeIndex < rg.vertices[curVertex].inEdges.size(); inEdgeIndex++) {
    inEdge = rg.vertices[curVertex].inEdges[inEdgeIndex];
    src = rg.edges[inEdge].src;
    if (rg.vertices[src].traversed == false) {
		  if (rg.edges[inEdge].traversed == false) {
  	    out << "\"" << rg.vertices[src].key << "\" -> \"" << rg.vertices[curVertex].key << "\" [label=\"" << rg.edges[inEdge].count << "\"];" << endl;
          rg.edges[inEdge].traversed = true;
			}
 		  if (rg.vertices[src].traversed == false) {
   	    RecursivePrintComponent(rg, src, out);
			}
    }
  }
  
}

int main(int argc, char* argv[]) {
	string rgInName, dotOutName;
	int minCoverage;

	if (argc < 3) {
		cout << "usage: removeTransitiveOverlaps in.rg out.dot [-separate]" << endl;
		exit(1);
	}

	rgInName    = argv[1];
	dotOutName   = argv[2];
	int argi = 3;
	bool printSeparate = false;
	while (argi < argc) {
		if (strcmp(argv[argi], "-separate") == 0) {
			printSeparate = true;
		}
		++argi;
	}
			

	RepeatGraph<string> rg;
	
	rg.ReadGraph(rgInName);

	VectorIndex vertexIndex;
	VectorIndex outEdgeIndex;
	VectorIndex edgeIndex;
	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++ ){
		rg.vertices[vertexIndex].traversed = false;
	}
	
	if (rg.edges.size() == 0) {
		cout << "LIKELY INVALID GRAPH. There are no edges." << endl;
		return 0;
	}
	//
	// At first, any edge that exists is connected to a vertex. This
	// will change as low coverage edges are deleted and replaced by
	// high coverage edges from the end of the array.
	//
	for (edgeIndex = 0; edgeIndex < rg.edges.size(); edgeIndex++) {
		rg.edges[edgeIndex].connected = true;
	}
  rg.Untraverse();
	if (printSeparate == false) {
		ofstream dotOut;
		CrucialOpen(dotOutName, dotOut, std::ios::out);
		
		dotOut << "digraph G{ " << endl;
		VectorIndex numEdges = rg.edges.size();
		for (edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
			dotOut << "\"" << rg.vertices[rg.edges[edgeIndex].src].key << "\" -> \"" << rg.vertices[rg.edges[edgeIndex].dest].key << "\"" 
						 << " [label=\"" << rg.edges[edgeIndex].count << "\"];" << endl;
		}
		dotOut<< "};" << endl;
	}
	else {
		bool foundUnprintedComponent = true;
		vertexIndex = 0;
		int componentIndex = 0;
		while (foundUnprintedComponent) {
			foundUnprintedComponent = false;
			for(; vertexIndex < rg.vertices.size(); vertexIndex++) {
				if (rg.vertices[vertexIndex].traversed == false) {
					//
					// found a new component.
					//
					cout << "printing comp: " << componentIndex << endl;
					stringstream componentNameStrm;
					componentNameStrm << dotOutName << "." << componentIndex << ".dot";
					string componentName = componentNameStrm.str();
					ofstream componentOut;
					CrucialOpen(componentName, componentOut, std::ios::out);
					componentOut << "digraph G { " << endl;
					RecursivePrintComponent(rg, vertexIndex, componentOut);			
				  componentOut << "}" << endl; 
				  componentOut.close();
	++componentIndex;
				}
			}
		}
	}
	return 0;
}
