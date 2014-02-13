#include "data/hdf/HDFCmpReader.h"
#include "data/hdf/HDFBasReader.h"
#include "CommandLineParser.h"
#include "FASTASequence.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "datastructures/alignment/Alignment.h"
#include "algorithms/alignment/AlignmentPrinter.h"
#include "algorithms/alignment/printers/XMLAlignmentPrinter.h"

int main(int argc, char* argv[]) {

	CommandLineParser clp;
	string cmpH5FileName, cmpXMLFileName;
	vector<int> holeNumbers;
	vector<string> patterns, refGroups;
  bool printAll = false;
	clp.RegisterStringOption("cmph5filename", &cmpH5FileName, "input cmp h5", false);
	clp.RegisterStringOption("cmpXMLfilename", &cmpXMLFileName, "output cmp xml", false);
	clp.RegisterPreviousFlagsAsHidden();
	clp.ParseCommandLine(argc, argv);

  HDFCmpReader<CmpAlignment> cmpReader;
  if (cmpReader.Initialize(cmpH5FileName) == 0) {
    cout << "Could not initialize the cmpH5 file, exiting." << endl;
    exit(0);
  }
  
  ofstream xmlOut;
  CrucialOpen(cmpXMLFileName, xmlOut, std::ios::out);
  
  unsigned int nAlignments = cmpReader.GetNAlignments();
  AlignmentCandidate<FASTASequence, FASTASequence> alignment;
  unsigned int i;

  for (i = 0; i < nAlignments; i++) {
    cmpReader.ReadAlignment(i, alignment);
    CompareXMLPrintAlignment(alignment, alignment.qAlignedSeq, alignment.tAlignedSeq, xmlOut);
  }

  xmlOut.close();


}

