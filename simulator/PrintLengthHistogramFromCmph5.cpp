#include "../common/data/hdf/HDFCmpReader.h"
#include "../common/datastructures/alignment/CmpFile.h"
#include "../common/utils.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: buildLengthHistogram.cpp input.cmp.h5 output.txt" << endl;
		exit(1);
	}

	string cmpH5FileName = argv[1];
	string outputFileName = argv[2];
	
	CmpFile cmpFile;
	HDFCmpReader<CmpAlignment> cmpFileReader;
	if (cmpFileReader.Initialize(cmpH5FileName) == 0) {
		cout << "ERROR, could not read the cmp file." << endl;
		exit(1);
	}
	ofstream outFile;
	CrucialOpen(outputFileName, outFile, std::ios::out);

	cmpFileReader.Read(cmpFile);
	vector<int> readLengths;
	int a;
	cout << "Processing " << cmpFile.alnInfo.alignments.size() << " alignments." << endl;
	for (a = 0; a < cmpFile.alnInfo.alignments.size(); a++ ) {
		readLengths.push_back(abs(cmpFile.alnInfo.alignments[a].GetQueryEnd() - cmpFile.alnInfo.alignments[a].GetQueryStart()));
	}
	sort(readLengths.begin(), readLengths.end());
	int r = 0;
	while(r < readLengths.size()) {
		int rn = r;
		while(rn < readLengths.size() and readLengths[rn] == readLengths[r]) { rn++;}
		outFile << readLengths[r] << " " << rn - r << endl;
		r = rn;
	}
}
