#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/utils.h"
#include "../common/NucConversion.h"

int main(int argc, char* argv[]) {

  string inFileName;
  string outFileName;
  ofstream outFile;
  inFileName = argv[1];
  outFileName = argv[2];

  CrucialOpen(outFileName, outFile, std::ios::out);

  FASTAReader reader;
  reader.Initialize(inFileName);

  FASTASequence forward;
  while(reader.GetNext(forward)) {
    int p;
    DNALength middle = forward.length / 2;
    for (p = 0; p < middle; p++) {
      Nucleotide right;
      DNALength rightIndex = forward.length - p - 1;
      right = forward.seq[rightIndex];
      forward.seq[rightIndex] = ReverseComplementNuc[forward.seq[p]];
      forward.seq[p] = ReverseComplementNuc[right];
    }
    if (middle % 2 != 0) {
      forward.seq[middle] = ReverseComplementNuc[forward.seq[middle]];
    }
    forward.PrintSeq(outFile);
  }
  
  return 0;
}
