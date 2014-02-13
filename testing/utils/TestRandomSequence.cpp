#include "algorithms/metagenomics/FindRandomSequence.h"
#include "FASTAReader.h"
#include "FASTASequence.h"

#include <string>

using namespace std;
int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: testRandomSequence genome.fa ntries " << endl;
    exit(0);
  }
  string inFile = argv[1];
  int nSamples = atoi(argv[2]);

  if (nSamples == 0) {
    return 0;
  }

  FASTAReader reader;
  reader.Initialize(inFile);
  vector<FASTASequence> genome;
  reader.ReadAllSequences(genome);
  
  int i;
  cout << "title pos" << endl;
  for (i = 0; i < nSamples; i++) {
    DNALength chrIndex, chrPos;
    FindRandomPos(genome, chrIndex, chrPos);
    cout << genome[chrIndex].title << " " << chrPos << endl;
  }

  return 0;
}
