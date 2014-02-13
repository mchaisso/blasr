#include "cdna/GencodeGFFFile.h"
#include "cdna/GencodeGFFGene.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "CommandLineParser.h"
#include <map>
#include <set>

using namespace std;
int main(int argc, char* argv[]) {
  string gencodeGffFileName, genesOutFileName;
  string geneType = "protein_coding";
  bool generateAll = false;
  if (argc < 3) {
    cout << "Usage: extractGenes gencodeGFFFile genesOutFileName [-geneType type [,alt_type] ]" << endl;
    exit(1);
  }

  gencodeGffFileName = argv[1];
  genesOutFileName   = argv[2];

  int argi = 4;
  string coordinatesFileName;
  set<string> geneTypes;

  while (argi < argc) {
    if (strcmp(argv[argi], "-geneType") == 0) {
      geneType = argv[++argi];
    }
    else {
      cout << "Bad option  " << argv[argi] << endl;
      exit(0);
    }
    ++argi;
  }

  vector<string> geneTypeList;
  ParseSeparatedList(geneType, geneTypeList, ',');
  int i;
  set<string> existingGenes;

  for (i = 0; i < geneTypeList.size(); i++) { geneTypes.insert(geneTypeList[i]); }

  ofstream outFile;
  CrucialOpen(genesOutFileName, outFile, std::ios::out);

  GencodeGFFFile gencodeFile;
  gencodeFile.ReadAll(gencodeGffFileName);
  
  vector<GencodeGFFGene> genes;
  IndexGencodeGenes(gencodeFile, genes, geneType);
  string curChromosome = "";

  for (i = 0; i < genes.size(); i++) {
    if (existingGenes.find(genes[i].geneName) != existingGenes.end()) {
      cout << "gene " << genes[i].geneName << " on chr " << genes[i].chromosome << endl;
    }
    existingGenes.insert(genes[i].geneName);
    genes[i].OrderExonsByStart();
  }

  for (i = 0; i < genes.size(); i++) {
    int e;
    if (curChromosome != genes[i].chromosome) {
      outFile << "START_CHROMOSOME" << endl;
      outFile << genes[i].chromosome << endl;
      curChromosome = genes[i].chromosome;
    }
    for (e = 0; e < genes[i].exons.size(); e++) {
      outFile << genes[i].geneName << " " << e 
              << " " << genes[i].exons[e]->start 
              << " " << genes[i].exons[e]->end 
              << " " << i << " NONE " << endl;
    }
  }
}
