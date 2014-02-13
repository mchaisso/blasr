#include "FASTAReader.h"
#include "FASTASequence.h"

#include <string>
#include <sstream>

int main(int argc, char* argv[]) {
  string inFileName, outFileName;
  int length;
  inFileName = argv[1];
  outFileName = argv[2];
  length = atoi(argv[3]);

  int argi = 4;
  int stride = 0;
  float coverage = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "-stride")) {
      stride = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-coverage")) {
      coverage = atof(argv[++argi]);
    }
    ++argi;
  }

  FASTAReader reader;
  reader.Initialize(inFileName);
  FASTASequence genome;
  reader.GetNext(genome);
  if (stride == 0 and coverage == 0) {
    cout << "error, must provide stride or coverage. " << endl;
    exit(0);
  }
  if (stride == 0) {
    stride = genome.length * coverage / length;
  }


  
  
