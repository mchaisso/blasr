#ifndef ASSEMBLY_PATH_OPERATIONS_H_
#define ASSEMBLY_PATH_OPERATIONS_H_

#include "RepeatGraph.h"
#include "Types.h"
void AdvanceSimplePath(RepeatGraph<string> &rg, UInt startVertex, UInt startEdge, UInt &endVertex, UInt &endEdge) {
	endVertex = rg.edges[startEdge].dest;
	endEdge   = startEdge;

	while(rg.vertices[endVertex].InDegree() == 1 and rg.vertices[endVertex].OutDegree() == 1) {
		endEdge   = rg.vertices[endVertex].outEdges[0];
		endVertex = rg.edges[endEdge].dest;
	}
}

void AdvanceAndStoreSimplePath(RepeatGraph<string> &rg, UInt startVertex, UInt startEdge, UInt &endVertex, UInt &endEdge, vector<UInt> &pathEdges) {
	endVertex = rg.edges[startEdge].dest;
	endEdge   = startEdge;
	pathEdges.push_back(endEdge);
	while(rg.vertices[endVertex].InDegree() == 1 and rg.vertices[endVertex].OutDegree() == 1) {
		endEdge   = rg.vertices[endVertex].outEdges[0];
		endVertex = rg.edges[endEdge].dest;
		pathEdges.push_back(endEdge);
	}
}


#endif
