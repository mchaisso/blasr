#include "FASTAReader.h"
#include "FASTASequence.h"
#include "algorithms/alignment/SWAlign.h"
#include <string>

int main(int argc, char* argv[]) {
	string vertexSequenceInFileName, vertexSequenceOutFileName;
	int lengthThreshold;
	float simThreshold;
	
	vertexSequenceInFileName = argv[1];
	lengthThreshold = atoi(argv[2]);
	simThreshold    = atof(argv[3]);
	vertexSequenceOutFileName = argv[4];

	

	

	
	

