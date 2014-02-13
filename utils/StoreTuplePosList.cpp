#include "tuples/DNATuple.h"
#include "tuples/TupleMetrics.h"
#include "tuples/DNATupleList.h"
#include "FASTAReader.h"
#include "FASTASequence.h"

#include <string>
#include <fstream>
#include <iostream>


using namespace std;
int main(int argc, char* argv[]) {

	string seqFileName;
	TupleMetrics tm;
	string outFileName;
	if (argc < 3) {
		cout << "usage: storeTuplePosList seqFile tupleSize outFile" << endl;
		return 0;
	}
	seqFileName = argv[1];
	tm.tupleSize = atoi(argv[2]);
	outFileName = argv[3];
	
	ofstream outFile;
	//	CrucialOpen(outFileName, outFile, std::ios::out| std::ios::binary);

	FASTAReader reader;
	reader.Init(seqFileName);
	FASTASequence seq;
	reader.GetNext(seq);
	//	vector<PositionDNATuple> 
	TupleList<PositionDNATuple>tuplePosList;
	tuplePosList.SetTupleMetrics(tm);
	//	StoreTuplePosList(seq, tm, tuplePosList);
	SequenceToTupleList(seq, tm, tuplePosList);
	tuplePosList.Sort();
	tuplePosList.WriteToFile(outFileName); //WriteTuplePosList(tuplePosList, tm.tupleSize, outFile);
	outFile.close();
	return 0;
}
