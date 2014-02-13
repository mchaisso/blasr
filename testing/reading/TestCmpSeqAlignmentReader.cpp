#include <string>
#include <iostream>
#include "datastructures/alignment/Alignment.h"
#include "algorithms/alignment/readers/CompareSequencesAlignmentReader.h"
#include "utils.h"

using namespace std;
int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "usage: testCmpSeqAlignmentReader cmpFile" << endl;
		exit(0);
	}
	string inFileName = argv[1];

	ifstream in;
	CrucialOpen(inFileName, in);
	CompSeqAlignment alignment;
	while(CompareSequencesAlignmentReader<CompSeqAlignment>::Read(in, alignment)) {
    int maxLength = alignment.qName.size();
    if (maxLength < alignment.tName.size()) {
      maxLength = alignment.tName.size();
    }
    string padding;
    int i;
    for (i = 0; i < maxLength; i++ ){ 
     padding.append(" ");
    }
    string qPadding, tPadding;
    for (i = alignment.qName.size(); i < maxLength; i++ ){
      qPadding.append(" ");
    }
    for (i = alignment.tName.size(); i < maxLength; i++ ){  
      tPadding.append(" ");
    }
   
 
		cout << alignment.qName << qPadding << " " << alignment.qString << endl;
    cout << padding << " " << alignment.alignString << endl;
		cout << alignment.tName << tPadding << " " << alignment.tString << endl;
	}
}
	
