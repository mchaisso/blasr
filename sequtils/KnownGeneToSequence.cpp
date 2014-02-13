#include "FASTASequence.h"
#include "FASTAReader.h"
#include "utils/StringUtils.h"
#include "utils.h"
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>

using namespace std;
int main(int argc, char* argv[]) {
  
  string genomeFileName, knownGenesFileName, knownGenesFastaFileName;
  if (argc < 4) {
    cout << "Usage: knownGeneToSequence genome.fasta knownGenes.txt output.fasta" << endl;
    exit(0);
  }
  genomeFileName = argv[1];
  knownGenesFileName = argv[2];
  knownGenesFastaFileName = argv[3];

  FASTAReader reader;
  reader.Initialize(genomeFileName);
  vector<FASTASequence> referenceSequences;
  reader.ReadAllSequences(referenceSequences);
  cout << "done reading sequences. " << endl;
  int i;
  map<string, int> titleToIndex;
  for (i = 0; i < referenceSequences.size(); i++) {
    titleToIndex[referenceSequences[i].title] = i;
  }

  ifstream knownGenesIn;
  CrucialOpen(knownGenesFileName, knownGenesIn, std::ios::in);

  ofstream rnaSeqOut;
  CrucialOpen(knownGenesFastaFileName, rnaSeqOut, std::ios::out);

  string line;
  while(getline(knownGenesIn, line)) {
    vector<unsigned int> exonStartPositions, exonEndPositions;
    stringstream lineStrm(line);
    // parse a line in the format:
    //uc009vit.3	chr1	-	14361	19759	14361	14361	9	14361,14969,15795,16606,16857,17232,17914,18267,18912,	14829,15038,15947,16765,17055,17742,18061,18366,19759,		uc009vit.3
    string title, chr;
    char strand;
    unsigned int start, end;
    unsigned int cdsStart, cdsEnd;
    int nExons;
    string exonStartString, exonEndString;
    lineStrm >> title >> chr >> strand >> start >> end >> cdsStart >> cdsEnd >> nExons >> exonStartString >> exonEndString;
    ParseSeparatedList(exonStartString, exonStartPositions);
    ParseSeparatedList(exonEndString, exonEndPositions);

    assert(exonStartPositions.size() == exonEndPositions.size());

    if (exonStartPositions.size() < 1) {
      cout <<"Error, no valid exon start positions for " << endl
           << line << endl;
      cout << "The exon string is " << exonStartString << endl;
      exit(0);
    }

    string rnaSeq;
    if (titleToIndex.find(chr) == titleToIndex.end()) {
      continue;
    }
    if (cdsStart == cdsEnd) { 
      continue;
    }
    /*
      cout << "ERROR. Could not find chromosome " << title << " in the genome." << endl;
      cout << "Valid chromosomes are: " << endl;
      map<string,int>::iterator mapIt;
      for (mapIt = titleToIndex.begin(); mapIt != titleToIndex.end(); ++mapIt) {
        cout << mapIt->first << endl;
      }
      return 1;
      }*/
    int chrIndex = titleToIndex[title];
    int e;
    for (e = 0; e < exonStartPositions.size() ; e++) {
      assert(exonEndPositions[e] > exonStartPositions[e]);
      rnaSeq.append((char*) &referenceSequences[chrIndex].seq[exonStartPositions[e]], exonEndPositions[e] - exonStartPositions[e]);
    }
    stringstream rnaSeqTitleStrm;
    rnaSeqTitleStrm << title << " " << " " << nExons << " " << exonStartString << " " << exonEndString;
    FASTASequence fastaRNASeq;
    
    fastaRNASeq.seq = (Nucleotide*) rnaSeq.c_str();
    fastaRNASeq.length = rnaSeq.size();
    if (strand == '+') {
      fastaRNASeq.CopyTitle(rnaSeqTitleStrm.str());
      fastaRNASeq.PrintSeq(rnaSeqOut);
    }
    else {
      FASTASequence rc;
      fastaRNASeq.MakeRC(rc);
      rc.CopyTitle(rnaSeqTitleStrm.str());
      rc.PrintSeq(rnaSeqOut);
    }
  }
}
