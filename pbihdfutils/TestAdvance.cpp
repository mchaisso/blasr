#include "SMRTSequence.h"
#include "data/hdf/HDFPlsReader.h"
#include "files/ReaderAgglomerate.h"
#include "utils/FileOfFileNames.h"
#include "utils.h"
#include <string>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {

	string plsFileName;
	int advance;

	if (argc <= 2) {
		cout << "usage: testAdvance file.pls.h5 advance " << endl;
		cout << "move 'advance' reads forward in a file." << endl;
		exit(1);
	}
	plsFileName = argv[1];
	advance = atoi(argv[2]);

	
	ReaderAgglomerate reader;
	reader.Initialize(plsFileName);
	
  SMRTSequence seq;
  int seqIndex = 0;
	int i;
	for (i = 0; i < 4; i++ ){
		seq.Free();
		reader.Advance(advance);
		reader.GetNext(seq);
	}
	seq.PrintSeq(cout);
}
