#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFBasWriter.h"
#include <string>
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
	
	string inFileName, outFileName;

	if (argc < 3) {
		cout << "usage: filterHDFPls in out idx1 [idx2 idx3]..." << endl;
		exit(1);
	}
	inFileName  = argv[1];
	outFileName = argv[2];
	
	vector<int> readIndices;
	int argi = 3;
	int minLength = 0;
	int minAvgQual = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minLength") == 0) {
			minLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minAvgQual") == 0) {
			minAvgQual = atoi(argv[++argi]);
		}
		++argi;
	}
	
	std::sort(readIndices.begin(), readIndices.end());
	HDFBasReader reader;
	HDFBasWriter writer;

	reader.Initialize(inFileName);
	writer.Initialize(outFileName, reader.GetMovieName(), reader.GetRunCode());
	
	int ri;
	int curReadIndex = 0;
	FASTQSequence seq;
	for (ri = 0; ri < readIndices.size(); ri++, curReadIndex++ ){
		reader.GetNext(seq);
		bool skipRead = false;
		if (seq.length < minLength) { skipRead = true;}
		if (seq.GetAverageQuality() < minAvgQual) { skipRead = true; }
		if (skipRead) { continue; }
		// all ok, write read out.
		writer.Write(seq);
	}


}
