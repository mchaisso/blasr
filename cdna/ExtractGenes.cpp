#include "cdna/GencodeGFFFile.h"
#include "cdna/GencodeGFFGene.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "CommandLineParser.h"
#include <map>


using namespace std;
int main(int argc, char* argv[]) {
  string gencodeGffFileName, genomeFileName, genesOutFileName;
  string geneType = "protein_coding";
  bool randomSplicing = false;
  int numRandomSplicing = 1;
  float pSkip = 0.5;
  if (argc < 4) {
    cout << "Usage: extractGenes gencodeGTFFile genomeFile genesOutFileName [-geneType type (protein_coding)] [-randomSplicing] [-numRandomSplicing n] [-pSkip prob (0-1, default:0.5)]" << endl;
    exit(1);
  }

  gencodeGffFileName = argv[1];
  genomeFileName     = argv[2];
  genesOutFileName   = argv[3];

  int argi = 4;
  string coordinatesFileName;

  while (argi < argc) {
    if (strcmp(argv[argi], "-geneType") == 0) {
      geneType = argv[++argi];
    }
    else if (strcmp(argv[argi], "-randomSplicing") == 0) {
      randomSplicing = true;
    }
    else if (strcmp(argv[argi], "-numRandomSplicing") == 0) {
      numRandomSplicing = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-pSkip") == 0) {
      pSkip = atof(argv[++argi]);
    }
    else {
      cout << "Bad option  " << argv[argi] << endl;
      exit(0);
    }
    ++argi;
  }

  coordinatesFileName = genesOutFileName;
  coordinatesFileName.append(".pos");
  FASTAReader reader;
  reader.Initialize(genomeFileName);

  ofstream outFile, coordsFile;
  CrucialOpen(genesOutFileName, outFile, std::ios::out);

  string coordsFileName = genesOutFileName + ".coords";
  CrucialOpen(coordsFileName, coordsFile, std::ios::out);

  vector<FASTASequence> referenceSequences;
  reader.ReadAllSequences(referenceSequences);
  int i;
  map<string, int> titleToIndex;
  for (i = 0; i < referenceSequences.size(); i++) {
    titleToIndex[referenceSequences[i].title] = i;
  }

  GencodeGFFFile gencodeFile;
  gencodeFile.ReadAll(gencodeGffFileName);
  
  vector<GencodeGFFGene> genes;
  IndexGencodeGenes(gencodeFile, genes, geneType);

  for (i = 0; i < genes.size(); i++) {
    genes[i].OrderExonsByStart();
  }

  int e;
  for (i = 0; i < genes.size(); i++) {
    FASTASequence geneSequence;
    geneSequence.CopyTitle(genes[i].geneName);
    if (titleToIndex.find(genes[i].chromosome) == titleToIndex.end()) {
      continue;
    }
    int chrIndex = titleToIndex[genes[i].chromosome];
    string sequence = "";
    //
    // Do nothing with 0 length exons.
    //
    if (genes[i].exons.size() == 0) {
      continue;
    }
    vector<FASTASequence> geneSequences;
    vector<GeneCoordinates> geneCoordinates;
    genes[i].GenerateGeneSequences(referenceSequences[chrIndex], geneSequences, geneCoordinates, randomSplicing);
    int gi;
    for (gi = 0; gi < geneSequences.size(); gi++) {
      if (genes[i].GetStrand() == '+') {
        geneSequences[gi].PrintSeq(outFile);
      }
      else {
        FASTASequence rc;
        geneSequences[gi].MakeRC(rc);
        rc.PrintSeq(outFile);
        rc.Free();
      }
      coordsFile << geneSequences[gi].title << " " << geneCoordinates[gi].chromosome << " " << geneCoordinates[gi].exonCoordinates.size() << " " << geneCoordinates[gi].strand;
      int i;
      for (i = 0; i < geneCoordinates[gi].exonCoordinates.size(); i++) {
        coordsFile << " " 
                   << geneCoordinates[gi].exonCoordinates[i].start << " "  
                   << geneCoordinates[gi].exonCoordinates[i].end << " ";
      }
      coordsFile << endl;
      geneSequences[gi].Free();
    }
    // 
    // No need to free the seq, since it is controlled by the string.
    //
  }
  coordsFile.close();
  
}
