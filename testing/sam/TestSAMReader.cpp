#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "datastructures/alignmentset/AlignmentSet.h"
#include "utils.h"
#include <string>
#include <iostream>

int main(int argc, char* argv[]) {
  
  if (argc < 2) {
    cout << "usage: testSamReader samFile" << endl;
    exit(1);
  }
  string samFileName = argv[1];
  
  AlignmentSet<> alignments;
  SAMReader<> samReader;
  samReader.Read(samFileName, alignments);
  
};
