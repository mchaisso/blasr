#include "data/hdf/HDFPlsReader.h"
#include "files/ReaderAgglomerate.h"
#include "utils/FileOfFileNames.h"
#include "SMRTSequence.h"
#include "utils.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

	string plsFileName, fastaOutName;

	if (argc < 2) {
		cout << "usage: pls2fasta  file.pls.h5 file.fasta " << endl;
		cout << "Print reads stored in hdf as fasta." << endl;
		exit(1);
	}
	vector<string> plsFileNames;
	plsFileName = argv[1];
	fastaOutName = argv[2];

	if (FileOfFileNames::IsFOFN(plsFileName)) {
		FileOfFileNames::FOFNToList(plsFileName, plsFileNames);
	}
	else {
		plsFileNames.push_back(plsFileName);
	}

	int plsFileIndex;
	for (plsFileIndex = 0; plsFileIndex < plsFileNames.size(); plsFileIndex++) {

		ReaderAgglomerate reader;
		reader.IgnoreCCS();
		reader.Initialize(plsFileNames[plsFileIndex]);

		ofstream fastaOut;
		CrucialOpen(fastaOutName, fastaOut);
	
		SMRTSequence seq;
		int seqIndex = 0;
		while (reader.GetNext(seq)) {
			seq.PrintQualSeq(fastaOut);
		}
	}
}
