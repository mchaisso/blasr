#ifndef ALGORITHMS_ALIGNMENT_READERS_ALIGNMENT_STATISTICS_READER_H_
#define ALGORITHMS_ALIGNMENT_READERS_ALIGNMENT_STATISTICS_READER_H_

#include <iostream>
#include <fstream>
using namespace std;
template<typename T_Alignment>
class AlignmentStatisticsReader {
 public:
	static int Read(ifstream &in, T_Alignment &alignment) {
		char qStrand;
		int qEnd, tEnd;
		char tStrand;
		if (!(in >> alignment.qName)) return 0;
		if (!(in >> alignment.qAlignLength )) return 0;
		if (!(in >> alignment.qPos )) return 0;
		if (!(in >> qEnd)) return 0;
		if (!(in >> qStrand )) return 0;
		if (qStrand == '-')
			alignment.qStrand = 1;
		else 
			alignment.qStrand = 0;
		if (!(in >> alignment.tName )) return 0;
		if (!(in >> alignment.tAlignLength )) return 0;
		if (!(in >> alignment.tPos )) return 0;
		if (!(in >> tEnd)) return 0;
		if (!(in >> tStrand )) return 0;
		if (tStrand == '-') 
			alignment.tStrand = 1;
		else
			alignment.tStrand = 0;
		if (!(in >> alignment.score )) return 0;
		if (!(in >> alignment.nMatch )) return 0;
		if (!(in >> alignment.nMismatch )) return 0;
		if (!(in >> alignment.nIns )) return 0;
		if (!(in >> alignment.nDel )) return 0;


		//
		// The alignments are output in the form (start end] rather than start, length
		// From the read above, alignments are stored in the [qt]Length field, so the 
		// length must be computed as an offset from the start ([qt]Pos).
		//
		alignment.qAlignLength -= alignment.qPos;
		alignment.tAlignLength -= alignment.tPos;
		return 1;
	}
};

#endif
