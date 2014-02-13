#ifndef ASSEMBLY_REPEAT_GRAPH_H_
#define ASSEMBLY_REPEAT_GRAPH_H_

#include <vector>
#include <string>
#include "RGVertex.h"
#include "RGEdge.h"
#include "utils.h"
#include "Types.h"
#include <assert.h>
int ReverseEdgeIndex(int index) { 
  return (2*(index/2) + !(index%2));
}

template<typename T_Key>
class RepeatGraph {
 public:
	vector<RGEdge>   edges;
	vector<RGVertex<T_Key> > vertices;
  void Untraverse() {
    UInt i;
    for (i = 0; i < edges.size(); i++ ){
    edges[i].traversed = false;
    }
    for (i = 0; i < vertices.size(); i++) {
    vertices[i].traversed = false;
    }
  }
 
	int Link(VectorIndex src, VectorIndex dest) {
		VectorIndex outEdgeIndex;
		VectorIndex srcOutIndex, destOutIndex;
		destOutIndex = srcOutIndex = 999999999;
		
		for (outEdgeIndex = 0; outEdgeIndex < vertices[src].outEdges.size(); outEdgeIndex++ ) {
			if (edges[vertices[src].outEdges[outEdgeIndex]].dest == dest) {
				srcOutIndex = outEdgeIndex;
				break;
			}
		}
		if (srcOutIndex >= vertices[src].outEdges.size()) {
			//
			// Create the new edge in the graph.
			//
			edges.push_back(RGEdge(src, dest));			
			//
			// Link the new edge into the graph.
			//
			vertices[src].outEdges.push_back(edges.size()-1);
			vertices[dest].inEdges.push_back(edges.size()-1);
			srcOutIndex = vertices[src].outEdges.size() - 1;
		}
		
		//
		// Update link counts for noise removal.
		//
		edges[vertices[src].outEdges[srcOutIndex]].count++;		
	}
	
	int LinkUndirected(VectorIndex src, VectorIndex dest) {
		//
		// An undirected link just uses the out edges.
		//
		VectorIndex outEdgeIndex;
		VectorIndex srcOutIndex, destOutIndex;
		destOutIndex = srcOutIndex = 999999999;
		
		for (outEdgeIndex = 0; outEdgeIndex < vertices[src].outEdges.size(); outEdgeIndex++ ) {
			if (edges[vertices[src].outEdges[outEdgeIndex]].dest == dest) {
				srcOutIndex = outEdgeIndex;
				break;
			}
		}
		if (srcOutIndex < vertices[src].outEdges.size()) {
			//
			// The two vertices are already linked.  Find that edge in the dest edge.
			//
			for (outEdgeIndex = 0; outEdgeIndex < vertices[dest].outEdges.size(); outEdgeIndex++ ) {
				if (edges[vertices[dest].outEdges[outEdgeIndex]].dest == src) {
					destOutIndex = outEdgeIndex;
					break;
				}
			}
		}
		if (srcOutIndex >= vertices[src].outEdges.size() or
				destOutIndex >= vertices[dest].outEdges.size()) {
			//
			// The two vertices are not yet linked.  Add the edge to source and dest.
			//
			edges.push_back(RGEdge(src, dest));
			vertices[src].outEdges.push_back(edges.size()-1);
			srcOutIndex = vertices[src].outEdges.size() - 1;
			
			edges.push_back(RGEdge(dest,src));
			vertices[dest].outEdges.push_back(edges.size()-1);
			destOutIndex = vertices[dest].outEdges.size() - 1;
		}
		
		//
		// Now update the weight of each link.
		//
		edges[vertices[src].outEdges[srcOutIndex]].count++;
		edges[vertices[dest].outEdges[destOutIndex]].count++;
	}
	
	void WriteGraph(ofstream &out) {
		out << vertices.size() << endl;
		VectorIndex vertexIndex;
		VectorIndex numVertices = vertices.size();
		for (vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
			out << vertices[vertexIndex];
		}
		out << edges.size() << endl;
		VectorIndex numEdges = edges.size();
		VectorIndex edgeIndex;
		for (edgeIndex = 0; edgeIndex < numEdges; edgeIndex++ ) {
			out << edgeIndex << " " << edges[edgeIndex];
		}
	}

	void WriteGraph(string &graphFileName) {
		ofstream out;
		CrucialOpen(graphFileName, out, std::ios::out);
		WriteGraph(out);
	}

	void ReadGraph(string &graphFileName) {
		ifstream in;
		CrucialOpen(graphFileName, in);

		VectorIndex numVertices;
		if (!(in >> numVertices)) {
			cout << "ERROR, could not read from graph" << endl;
			exit(1);
		}
		VectorIndex vertexIndex;
		vertices.resize(numVertices);
		for (vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
			in >> vertices[vertexIndex];
		}

		VectorIndex numEdges;
		in >> numEdges;
		edges.resize(numEdges);
		VectorIndex edgeIndex;
		VectorIndex tmpEdgeIndex;
		for (edgeIndex = 0; edgeIndex < numEdges; edgeIndex++ ) {
			in >> tmpEdgeIndex >> edges[edgeIndex];
		}
	}
	
};


#endif
