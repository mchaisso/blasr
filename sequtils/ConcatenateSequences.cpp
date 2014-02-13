#include <vector>
#include <string>
#include "FASTASequence.h"
#include "FASTAReader.h"
#include "NucConversion.h"
#include "utils.h"
#include "CommandLineParser.h"
#include "datastructures/metagenome/SequenceIndexDatabase.h"
#include "utils/FileUtils.h"

int main(int argc, char* argv[]) {
  CommandLineParser clp;
  string fastaOutName;
	vector<string> inFiles;  
  string indexFileName;
  clp.RegisterStringOption("outFile", &fastaOutName, "Write to this file.", true);
  clp.RegisterStringListOption("inFiles", &inFiles, "Read from these files.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterStringOption("writeIndex", &indexFileName, "Print a map of read name to start,end coordinates of the output file.");

  clp.ParseCommandLine(argc, argv);

  if (FileExists(fastaOutName)) {
    cout << "ERROR. The file " << fastaOutName << " already exists.  This will not be overwritten in case there is an error and you lose precious data."<< endl;
    exit(1);
  }

  ofstream indexFileOut;
  if (indexFileName!="") {
    CrucialOpen(indexFileName, indexFileOut, std::ios::out);
  }

  SequenceIndexDatabase<FASTASequence> seqDB;
  
	//
	// Read the suffix array to modify.
	//
	
	ofstream fastaOut;
	CrucialOpen(fastaOutName, fastaOut);
	VectorIndex inFileIndex;
	int seqIndex = 0;
	FASTASequence seq;
	for (inFileIndex = 0; inFileIndex < inFiles.size(); ++inFileIndex) {
		FASTAReader reader;
		reader.Init(inFiles[inFileIndex]);
    reader.ReadAllSequencesIntoOne(seq, &seqDB);
    if (inFileIndex == 0) {
      seq.PrintSeq(fastaOut);
    }
    else {
      // Print without the title
      ((DNASequence)seq).PrintSeq(fastaOut);
    }
	}
  
  if (indexFileName != "") {
    int i;
    if (seqDB.growableName.size() > 1) {
      for (i = 0; i < seqDB.growableName.size(); i++) {
        indexFileOut << seqDB.growableName[i] << " " << seqDB.growableSeqStartPos[i] << " " << seqDB.growableSeqStartPos[i+1] << endl;
      }
    }
  }
	return 0;
}
