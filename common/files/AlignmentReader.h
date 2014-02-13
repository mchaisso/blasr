#ifndef FILES_ALIGNMENT_READER_H_
#define FILES_ALIGNMENT_READER_H_

#include "../algorithms/alignment/AlignmentFormats.h"
#include "../FASTASequence.h"
#include "../utils.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

template<typename T_Alignment>
class AlignmentReader {
	public:
	ifstream in;
	AlignmentPrintFormat fileFormat;

	int Init(string alignmentFileName) {
		CrucialOpen(alignmentFileName, in);
		fileFormat = DetermineAlignmentFormatType(alignmentFileName);
	}

	AlignmentPrintFormat DetermineAlignmentFormatType(string alignmentFileName) {
		string suffix;
		string::size_type dotPos = alignmentFileName.find_last_of('.');
		if (dotPos == alignmentFileName.npos) {
			return NOFORMAT;
		}
		else {
			suffix = alignmentFileName.substr(dotPos+1);
			if (suffix == "rm0" or suffix == "m0") {
				return StickPrint;
			}
			else if (suffix == "rm1" or suffix == "m1") {
				return SummaryPrint;
			}
			else if (suffix == "rm4" or suffix == "m4") {
				return Interval;
			}
			else if (suffix == "rm5" or suffix == "m5") {
				return CompareSequencesParsable;
			}
			else {
				return NOFORMAT;
			}
		}
	}

	int GetNext(T_Alignment &alignment) {
		switch(fileFormat){ 
		case(StickPrint):
			cout << "ERROR, reading of stick alignments is not implemented NOR recommended!!" << endl;
			exit(1);
		case(SummaryPrint):
			return ReadSummaryAlignment(alignment);
			break;  // add break so that if return is removed, bugs dont arise
		case(Interval):
			return ReadIntervalAlignment(alignment);
			break; 
		case(CompareSequencesParsable):
			return ReadCompareSequencesParsableAlignment(alignment);
			break;
		case(NOFORMAT):
			return 0; // this shouldn't happen maybe assert?
			break;
		}
	}

	int ReadStickPrintAlignment(T_Alignment &alignment) {
		return 0;
	}
	
	int ReadSummaryAlignment(T_Alignment &alignment) {
		cout << "ReadSummaryAlignment not yet implemented." << endl;
		return 0;
	}

	int ReadIntervalAlignment(T_Alignment &alignment) {
		//
		// Read in values one by one to make sure they all exist.
		//
		if (! (in >> alignment.qName)) return 0;
		if (! (in >> alignment.tName)) return 0;
		if (! (in >> alignment.score)) return 0;
		if (! (in >> alignment.probScore)) return 0;
		if (! (in >> alignment.sumQVScore)) return 0;
		if (! (in >> alignment.pctSimilarity)) return 0;
		if (! (in >> alignment.qStrand)) return 0;
		if (! (in >> alignment.qPos)) return 0;
		DNALength qEnd, tEnd;
		if (! (in >> qEnd)) return 0;
		if (! (in >> alignment.qLength)) return 0;
		if (! (in >> alignment.tStrand)) return 0;
		if (! (in >> alignment.tPos)) return 0;
		if (! (in >> tEnd)) return 0;
		if (! (in >> alignment.tLength)) return 0;

		//
		// Assign length considering possible swap of directions of strand.
		//			
		return 1;
	}
		
	int ReadCompareSequencesParsableAlignment(T_Alignment &alignment) {
		cout << "ReadCompareSequencesParsable is not yet implemented." << endl;
	}

	void ReadAllAlignments(vector<T_Alignment> &alignments) { 
		T_Alignment alignment;
		while (ReadIntervalAlignment(alignment)) {
			alignments.push_back(alignment);
		}
	}
};



#endif
