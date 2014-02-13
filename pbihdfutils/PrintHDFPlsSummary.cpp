#include "data/hdf/HDFBasReader.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

	string basFileName;

	if (argc < 2) {
		cout << "usage: bas2  file.bas.h5 " << endl;
		cout << "bas2 = bas ls" << endl;
		cout << " -t --titleindices  Print the titles of reads and their indices in the file." << endl;
		exit(1);
	}
	basFileName = argv[1];
	int argi = 1;
	bool printSummary      = true;
	bool printTitleIndices = false;
	while (argi < argc) {
		if (strcmp(argv[argi], "-t") == 0 or
				strcmp(argv[argi], "--titleindices") == 0) {
			printTitleIndices = true;
			printSummary = false;
		}
		++argi;
	}
	HDFBasReader reader;
	reader.Initialize(basFileName);
	int numReads = reader.GetNumReads();
	cout << "Num reads: " << numReads << endl;
	vector<int> readLengths;
	reader.GetAllReadLengths(readLengths);
	DNALength totalLength = 0;
	VectorIndex i;
	int numAbove100 = 0;
	for (i = 0; i< readLengths.size(); i++ ){
		totalLength += readLengths[i];
		if (readLengths[i] > 100) {
			++numAbove100;
		}
	}

	if (printSummary) {
		if (numReads > 0 ){
			cout << "Average read length: " << totalLength / (1.0*numReads) << endl;
		}
		cout << "A total of " << numAbove100  << " have length > 100" << endl;

		return 0;
	}
	if (printTitleIndices) {
		FASTQSequence seq;
		int seqIndex = 0;
		while (reader.GetNext(seq)) {
			cout << seqIndex << " " << seq.title << " " << seq.length << endl;
			++seqIndex;
		}
	}
}
