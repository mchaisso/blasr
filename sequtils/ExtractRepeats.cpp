#include "../common/datastructures/repmask/RepmaskTable.h"
#include "../common/utils.h"
#include "../common/FASTAReader.h"
#include <iostream>
#include <fstream>
#include <sstream>


void PrintUsage() {
  cout << "usage: extractrep  genome.fa genome.fa.out repatOut.fasta [-repeat rep] [-family fam]" << endl;
}


int main(int argc, char* argv[]) {
  string genomeFileName, dotoutFileName, repeatOutFileName, repeatName, repeatFamily;
  if (argc < 4) {
    PrintUsage();
    exit(1);
  }
  genomeFileName = argv[1];
  dotoutFileName = argv[2];
  repeatOutFileName = argv[3];

  int argi = 4;
  while (argi < argc) {
    if (strcmp(argv[argi], "-name") == 0) {
      repeatName = argv[++argi];
    }
    else if (strcmp(argv[argi], "-family") == 0) {
      repeatFamily = argv[++argi];
    }
    else {
      cout << "ERROR! Option " << argv[argi] << " is not valid." << endl;
      PrintUsage();
      exit(1);
    }
    ++argi;
  }

  FASTAReader reader;
  reader.Initialize(genomeFileName);
  FASTASequence genome;
  reader.GetNext(genome);

  RepmaskTable repeatTable;
  repeatTable.Read(dotoutFileName);

  ofstream repeatSeqOut;
  CrucialOpen(repeatOutFileName, repeatSeqOut, std::ios::out);
  
  int i;
  for (i = 0; i < repeatTable.size(); i++) {
    bool printRepeat = false;
    if (repeatName != "" and repeatName == repeatTable[i].matchingRepeat ) {
      printRepeat = true;
    }
    else if (repeatFamily != "" and repeatFamily == repeatTable[i].repeatClass) {
      printRepeat = true;
    }
    if (printRepeat) {
      FASTASequence seq;
      stringstream sstrm;
      sstrm << repeatTable[i].matchingRepeat << "/" << repeatTable[i].repeatClass 
            << "/" << repeatTable[i].genomeStartPos << "/" << repeatTable[i].genomeEndPos;
      seq.CopyTitle(sstrm.str());
      seq.seq = &genome.seq[repeatTable[i].genomeStartPos-1];
      seq.length = repeatTable[i].genomeEndPos - repeatTable[i].genomeStartPos;
      seq.PrintSeq(repeatSeqOut);
    }
  }
}

