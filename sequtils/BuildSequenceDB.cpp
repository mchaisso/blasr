#include <string>
#include "../common/datastructures/metagenome/SequenceIndexDatabase.h"
#include "../common/CommandLineParser.h"
#include "../common/FASTAReader.h"
#include "../common/utils/FileOfFileNames.h"
#include "../common/utils.h"

using namespace std;

int main(int argc, char* argv[]) {

	CommandLineParser clp;
	string fastaFileName, indexFileName;
	vector<string> fastaFileNames;
	vector<string> opts;
	clp.SetProgramName("bsdb");
	clp.SetProgramSummary("Build an index database on a file of sequences.\n"
												" The index is used to map to reads given alignment positions.\n");
	clp.RegisterStringOption("fasta", &fastaFileName, "A file with sequences to build an index.");
	clp.RegisterStringOption("index", &indexFileName, "The index file.");
	clp.RegisterPreviousFlagsAsHidden();

	clp.ParseCommandLine(argc, argv, opts);

	ifstream fastaIn;
	ofstream indexOut;

	if (FileOfFileNames::IsFOFN(fastaFileName)) {
		FileOfFileNames::FOFNToList(fastaFileName, fastaFileNames);
	}
	else {
		fastaFileNames.push_back(fastaFileName);
	}

	CrucialOpen(indexFileName, indexOut, std::ios::out | std::ios::binary);
	SequenceIndexDatabase<FASTASequence> seqDB;
		
	int fileNameIndex;
	for (fileNameIndex = 0; fileNameIndex < fastaFileNames.size(); fileNameIndex++){ 
		FASTAReader reader;
		FASTASequence seq;
		reader.Init(fastaFileNames[fileNameIndex]);
		int i = 0;
		while (reader.GetNext(seq)) {
			seqDB.AddSequence(seq);
			i++;
		}
	}
	seqDB.Finalize();
	seqDB.WriteDatabase(indexOut);
	return 0;
}
