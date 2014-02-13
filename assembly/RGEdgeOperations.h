#ifndef REPEAT_GRAPH_EDGE_OPERATIONS_H_
#define REPEAT_GRAPH_EDGE_OPERATIONS_H_

#include "RepeatGraph.h"
#include <set>
#include <vector>
using namespace std;

template<typename T_Vertex, typename T_Edge>
void UnlinkDirected(vector<T_Vertex> &vertices, vector<T_Edge> &edges, VectorIndex edgeIndex) {
	VectorIndex src, dest;
	src  = edges[edgeIndex].src;
  dest = edges[edgeIndex].dest;
	VectorIndex outEdgeIndex, inEdgeIndex;

	//
	// Remove the outEdge from src
	//
	for (outEdgeIndex = 0; outEdgeIndex < vertices[src].outEdges.size(); outEdgeIndex++) {
		if (vertices[src].outEdges[outEdgeIndex] == edgeIndex) {
			break;
		}
	}
	// This edge must exist in the out edge list of this vertex.
	assert(outEdgeIndex < vertices[src].outEdges.size());

	for(; outEdgeIndex < vertices[src].outEdges.size()-1; outEdgeIndex++ ){
		vertices[src].outEdges[outEdgeIndex] = vertices[src].outEdges[outEdgeIndex+1];
	}
	vertices[src].outEdges.resize(vertices[src].outEdges.size()-1);


	//
	// Remove the in edge from dest
	//
	for (inEdgeIndex = 0; inEdgeIndex < vertices[dest].inEdges.size(); inEdgeIndex++) {
		if (vertices[dest].inEdges[inEdgeIndex] == edgeIndex) {
			break;
		}
	}
	// This edge must exist in the inEdges list of this vertex.
	assert(inEdgeIndex < vertices[dest].inEdges.size());

	for(; inEdgeIndex < vertices[dest].inEdges.size()-1; inEdgeIndex++ ){
		vertices[dest].inEdges[inEdgeIndex] = vertices[dest].inEdges[inEdgeIndex+1];
	}
	vertices[dest].inEdges.resize(vertices[dest].inEdges.size()-1);
	edges[edgeIndex].connected = false;
}


//
// Undirected version of above.
//
template<typename T_Vertex, typename T_Edge>
void UnlinkSrc(vector<T_Vertex> &vertices, vector<T_Edge> &edges, VectorIndex edgeIndex) {
	VectorIndex src;
	src = edges[edgeIndex].src;
	VectorIndex outEdgeIndex;
	for (outEdgeIndex = 0; outEdgeIndex < vertices[src].outEdges.size(); outEdgeIndex++) {
		if (vertices[src].outEdges[outEdgeIndex] == edgeIndex) {
			break;
		}
	}
	// This edge must exist in the out edge list of this vertex.
	assert(outEdgeIndex < vertices[src].outEdges.size());

	for(; outEdgeIndex < vertices[src].outEdges.size()-1; outEdgeIndex++ ){
		vertices[src].outEdges[outEdgeIndex] = vertices[src].outEdges[outEdgeIndex+1];
	}
	vertices[src].outEdges.resize(vertices[src].outEdges.size()-1);
	edges[edgeIndex].connected = false;
}

template<typename T_Vertex, typename T_Edge>
void UpdateOutEdge(vector<T_Vertex> &vertices, vector<T_Edge> &edges,VectorIndex prevEdgeIndex, VectorIndex newEdgeIndex) {
	VectorIndex src;
	// The new edge must be already updated.
	src = edges[newEdgeIndex].src;
	
	VectorIndex outEdgeIndex;
	for (outEdgeIndex = 0; outEdgeIndex < vertices[src].outEdges.size(); outEdgeIndex++) {
		if (vertices[src].outEdges[outEdgeIndex] == prevEdgeIndex) {
			break;
		}
	}
	// This edge must exist in the out edge list of this vertex.
	assert(outEdgeIndex < vertices[src].outEdges.size());

	// Relink the old out edge slot to the new one.
	vertices[src].outEdges[outEdgeIndex] = newEdgeIndex;
}

template<typename T_Vertex, typename T_Edge>
void UpdateInEdge(vector<T_Vertex> &vertices, vector<T_Edge> &edges,VectorIndex prevEdgeIndex, VectorIndex newEdgeIndex) {
	VectorIndex dest;
	// The new edge must be already updated.
	dest = edges[newEdgeIndex].dest;
	
	VectorIndex inEdgeIndex;
	for (inEdgeIndex = 0; inEdgeIndex < vertices[dest].inEdges.size(); inEdgeIndex++) {
		if (vertices[dest].inEdges[inEdgeIndex] == prevEdgeIndex) {
			break;
		}
	}
	// This edge must exist in the in edge list of this vertex.
	assert(inEdgeIndex < vertices[dest].inEdges.size());

	// Relink the old in edge slot to the new one.
	vertices[dest].inEdges[inEdgeIndex] = newEdgeIndex;
}

#endif

template<typename T_Vertex, typename T_Edge>
void RemoveUnconnectedEdges(vector<T_Vertex> &vertices, vector<T_Edge> &edges) {
	VectorIndex numEdges = edges.size();
	VectorIndex edgeIndex;
	for (edgeIndex = 0; edgeIndex < numEdges; ) {
		if (edges[edgeIndex].connected == true) {
			// This edge is fine.
			edgeIndex++;
		}
		else {
			// This edge needs to be deleted. Find the first edge at the end that
			// will not be deleted anyway, and move it here.

			while (numEdges > edgeIndex and edges[numEdges-1].connected == false) {
				UnlinkDirected(vertices, edges, numEdges-1);
				numEdges--;
			}

			//
			// If exhausted all edges, just break since all are deleted.
			//
			if (numEdges == edgeIndex) {
				continue;
			}
			VectorIndex src  = edges[edgeIndex].src;
			VectorIndex dest = edges[edgeIndex].dest;
			//
			// Get rid of this low coverage edge.
			UnlinkDirected(vertices, edges, edgeIndex);

			//
			// Pack in one from a higher coverage.
			edges[edgeIndex] = edges[numEdges-1];

			//
			// Update the connecting vertex.
			UpdateOutEdge(vertices, edges, numEdges-1, edgeIndex);
			UpdateInEdge(vertices, edges, numEdges-1, edgeIndex);			
			--numEdges;

			++edgeIndex;
		}
	}
	edges.resize(numEdges);
}

template<typename T_Vertex, typename T_Edge>
void MarkEdgePairsForRemoval(set<pair<VectorIndex, VectorIndex> > &srcDestToRemove, 
												 vector<T_Vertex> &vertices, 
												 vector<T_Edge> &edges) {

  set<pair<VectorIndex,VectorIndex> >::iterator etrFirst, etrLast;
  etrFirst = srcDestToRemove.begin();
	etrLast  = srcDestToRemove.end();

	//
	// Mark edges that should be removed.
	//
	VectorIndex edgeIndex, outEdgeIndex;
  while(etrFirst != etrLast) {
		VectorIndex src, dest;
		src = etrFirst->first;
		dest = etrFirst->second;
		for (outEdgeIndex = 0; outEdgeIndex < vertices[src].outEdges.size(); outEdgeIndex++) {
			edgeIndex = vertices[src].outEdges[outEdgeIndex];
			if ( edges[edgeIndex].dest == dest) {
				edges[edgeIndex].connected = false;
				break;
			}
		}
		++etrFirst;
  }
}
