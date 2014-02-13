#include "data/hdf/HDFCmpFile.h"
#include "data/hdf/HDFBasReader.h"
#include "CommandLineParser.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "datastructures/alignment/Alignment.h"
#include "algorithms/alignment/AlignmentPrinter.h"

		
		
		

int main(int argc, char* argv[]) {


	CommandLineParser clp;
	string cmpFileName;
	vector<int> holeNumbers;
	vector<string> patterns, refGroups;
  bool printAll = false;
	clp.RegisterStringOption("cmph5filename", &cmpFileName, "input cmp h5", false);
	clp.RegisterPreviousFlagsAsHidden();
	clp.RegisterIntListOption("holeNumbers", &holeNumbers, "hole numbers to print alignments", false);
	clp.RegisterStringListOption("pattern", &patterns, "patterns to search read names to print alignments", false);	
  clp.RegisterFlagOption("all", &printAll, "Just print all alignments.", false);
  clp.RegisterStringListOption("refgroups", &refGroups, "Reference groups to print.", false);
	clp.ParseCommandLine(argc, argv);

	
	CmpFile cmpFile;
	
	/*
	 * These readers pull information from the same pls file.
	 */
	HDFCmpFile<CmpAlignment> hdfcmpFile;

	if (hdfcmpFile.Initialize(cmpFileName) == 0) {
		cout << "ERROR, could not open the cmp file." << endl;
		exit(1);
	}
	
	hdfcmpFile.Read(cmpFile);
	
	int alignmentIndex;
	for (alignmentIndex = 0; alignmentIndex < cmpFile.alnInfo.alignments.size(); alignmentIndex++) {
		int alnHoleNumber;
		alnHoleNumber = cmpFile.alnInfo.alignments[alignmentIndex].GetHoleNumber();
		int hi;
    bool printThisAlignment = false;

    //
    // Read the alignment string.  All alignments 
    //
    int refGroupId = cmpFile.alnInfo.alignments[alignmentIndex].GetRefGroupId();
    int alnGroupId  = cmpFile.alnInfo.alignments[alignmentIndex].GetAlnGroupId();

    int refGroupIndex = hdfcmpFile.refGroupIdToArrayIndex[refGroupId];
    string readGroupName = hdfcmpFile.alnGroupIdToReadGroupName[alnGroupId];
    int readGroupIndex = hdfcmpFile.refAlignGroups[refGroupIndex]->experimentNameToIndex[readGroupName];

    string refGroupPath = cmpFile.refGroup.path[refGroupIndex];

		for (hi = 0; hi < holeNumbers.size(); hi++) {
			if (alnHoleNumber == holeNumbers[hi]) {
        printThisAlignment = true;
        break;
      }
    }
    int ri;
    for (ri = 0; ri < refGroups.size(); ri++) {
      if (refGroups[ri] == refGroupPath) {
        printThisAlignment = true;
        break;
      }
    }


    if (printThisAlignment or printAll) {
      unsigned int alignStartIndex, alignEndIndex;
      UInt offsetBegin, offsetEnd;
		
      string   refSequence;
      string   readSequence;
      vector<unsigned char> byteAlignment;

      offsetBegin = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
      offsetEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd();
      int alignedSequenceLength = offsetEnd - offsetBegin;
      if (alignedSequenceLength >= 0) {
        refSequence.resize(alignedSequenceLength);
        byteAlignment.resize(alignedSequenceLength);
      }
	
      
      hdfcmpFile.refAlignGroups[refGroupIndex]->readGroups[readGroupIndex]->alignmentArray.Read(offsetBegin, 
                                                                                               offsetEnd, 
                                                                                               &byteAlignment[0]);

      readSequence.resize(byteAlignment.size());
      refSequence.resize(byteAlignment.size());

      ByteAlignmentToQueryString(&byteAlignment[0], byteAlignment.size(), &readSequence[0]);
      ByteAlignmentToRefString(&byteAlignment[0], byteAlignment.size(), &refSequence[0]);				
      string ungappedRead, ungappedRef;
      RemoveGaps(readSequence, ungappedRead);
      RemoveGaps(refSequence, ungappedRef);
      Alignment alignment;
      GappedStringsToAlignment(readSequence, refSequence, alignment);
      DNASequence qAlignedSeq, rAlignedSeq;
      qAlignedSeq.seq = (Nucleotide*) &ungappedRead[0];
      qAlignedSeq.length = ungappedRead.size();
      rAlignedSeq.seq = (Nucleotide*) &ungappedRef[0];
      rAlignedSeq.length = ungappedRef.size();
				
      int qStart = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart();
      int tStart = cmpFile.alnInfo.alignments[alignmentIndex].GetRefStart();
      stringstream sstrm;
      sstrm << alnHoleNumber << "/" << qStart << "_" << cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd();
      alignment.qName = sstrm.str();
      StickPrintAlignment(alignment, qAlignedSeq, rAlignedSeq, cout, qStart, tStart);
				
    }
  }
}

