#include "cdna/GeneDB.h"
#include "cdna/GencodeGFFFile.h"
#include "cdna/GencodeGFFGene.h"
#include "utils.h"
#include <map>
#include <set>


int main(int argc, char* argv[]) {

  string alignmentFileName, geneDBFileName;

  if (argc < 3) {
    cout << "usage: printDuplicates alignemntFile genedb " << endl;
    exit(0);
  }

  alignmentFileName = argv[1];
  geneDBFileName    = argv[2];

  GeneDB genedb;
  genedb.Read(geneDBFileName);
  genedb.IndexChromosomes();

  ifstream alignmentIn;
  CrucialOpen(alignmentFileName, alignmentIn, std::ios::in);
  
  while(alignmentIn) {
    string line;
    getline(alignmentIn, line);
    if (line == "###") {
      // found the end of an entry
      

    }
    else {
      string chrName, genome, source;
      int start, end, identity;
      char a, strand, b;
      string annotationString;
      stringstream strm(line);
      strm >> chrName >> genome >> source >> start >> end >> a >> strand >> b >> annotationString;
      vector<string> annotations;
      if (source != "exon") {
        continue;
      }
      ParseSeparatedList(annotationString, annotations, ';');

      GeneDBChromosome *chromosomePtr;
      chromosomePtr = genedb.Find(chrName);
      if (chromosomePtr == NULL) {
        cout << "chromosome " << chromosomePtr << " not found in database." << endl;
        continue;
      }
      int exonIndex;
      if ( chromosomePtr->LookupIndexByStart(start, exonIndex) ) {
        if (chromosomePtr->exons[exonIndex].start <= start and 
            chromosomePtr->exons[exonIndex].end   >= end) {
          chromosomePtr->exons[exonIndex].Print(cout);
        }
      }
      else {
        
      }
    }
  }
}
  

