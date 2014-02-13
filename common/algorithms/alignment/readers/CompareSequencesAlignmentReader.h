#ifndef ALGORITHMS_ALIGNMENT_READERS_COMPARE_SEQUENCES_ALIGNMENT_READER_H_
#define ALGORITHMS_ALIGNMENT_READERS_COMPARE_SEQUENCES_ALIGNMENT_READER_H_

#include "AlignmentStatisticsReader.h"
using namespace std;
template<typename T_Alignment>
class CompareSequencesAlignmentReader {
 public:
	static int Read(ifstream &in, T_Alignment &alignment) {
    AlignmentStatisticsReader<T_Alignment>::Read(in, alignment);
		if (!(in >> alignment.qString)) return 0;
		if (!(in >> alignment.alignString)) return 0;
		if (!(in >> alignment.tString)) return 0;
		return 1;
	}

	
};


#endif
