#ifndef DATASTRUCTURES_METAGENOME_FASTA_TITLE_DICTIONARY_H_
#define DATASTRUCTURES_METAGENOME_FASTA_TITLE_DICTIONARY_H_

#include <string>
#include <map>
#include "../../utils.h"
#include "../../FASTASequence.h"
#include "../../Types.h"
using namespace std;

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

typedef map<string, unsigned int> TitleMap;
class FASTATitleDictionary {
	
 public:
	TitleMap titleMap;
	int curSeq;

	FASTATitleDictionary() { 
		curSeq = 0;
	}

	void AddAllSequences(vector<FASTASequence> &sequences){
		VectorIndex seqIndex;
		for (seqIndex = 0; seqIndex < sequences.size(); seqIndex++) {
			AddSequence(sequences[seqIndex]);
		}
	}

	void AddSequence(FASTASequence &seq) {
		string seqName = seq.GetName();
		titleMap[seqName] = curSeq; 
		curSeq++;
	}

	int LookupSequence(string seqName, int &seqIndex) {
		TitleMap::iterator it;
		it = titleMap.find(seqName);
		if ( it != titleMap.end()) {
			seqIndex = it->second;
			return 1;
		}
		else {
			seqIndex = 0;
			return 0;
		}
	}
};


#endif
