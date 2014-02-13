#include <string>
#include "FASTASequence.h"
#include "datastructures/metagenome/SequenceIndexDatabase.h"
#include "CommandLineParser.h"
#include "utils.h"

int main(int argc, char* argv[]) {
	CommandLineParser clp;
	
	string indexDBName;
	bool   printIndex = false;
	int    searchIndex;

	//
	// Configure the command line.
	//
	clp.SetProgramName("testseqdb");
	clp.SetProgramSummary("test the sequence db.\n");
	clp.RegisterStringOption("indexdb", &indexDBName, "The index to test.");
	clp.RegisterPreviousFlagsAsHidden();
	clp.RegisterFlagOption("p", &printIndex, "Print the start position of each read.");
	clp.RegisterIntOption("i", &searchIndex, "The index to search for", CommandLineParser::NonNegativeInteger, true);

	clp.ParseCommandLine(argc, argv);

	SequenceIndexDatabase<FASTASequence> seqDB;
	ifstream in;
	CrucialOpen(indexDBName, in, std::ios::in | std::ios::binary);
	seqDB.ReadDatabase(in);

	if (printIndex) {
		int i;
		for (i = 0; i < seqDB.nSeqPos - 1; i++) { 
			cout << i << " " << seqDB.seqStartPos[i+1] << " " << seqDB.names[i] << endl;
		}
	}

	int dbPos = seqDB.SearchForIndex(searchIndex);
	if (dbPos >= 0) {
		cout << "searchIndex: " << searchIndex << " " << dbPos << " " << seqDB.seqStartPos[dbPos] << " " << seqDB.names[dbPos-1] << endl;
	}
};
	
