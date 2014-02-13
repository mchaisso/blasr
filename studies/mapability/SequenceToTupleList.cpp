#include "tuples/DNATupleList.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <string>
#include "tuples/TupleMetrics.h"
using namespace std;

int main(int argc, char* argv[]) {

	string inFileName, outFileName;
	TupleMetrics tm;
	inFileName = argv[1];
	tm.tupleSize = atoi(argv[2]);
	outFileName =argv[3];

	ofstream out;
	CrucialOpen(outFileName, out, std::ios::out|std::ios::binary);

	DNASequence in;
	FASTAReader reader;
	reader.Init(inFileName);
	FASTASequence seq;
	reader.ReadAllSequencesIntoOne(seq);
	TupleList<DNATuple> tupleList;
	SequenceToTupleList(seq, tm, tupleList);
	out.write((char*) &tupleList[0], sizeof(DNATuple) * tupleList.size());
	return 0;

}
