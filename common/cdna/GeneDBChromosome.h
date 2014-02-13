#ifndef CDNA_GENE_DB_CHROMOSOME_H_
#define CDNA_GENE_DB_CHROMOSOME_H_
#include "GeneDBExon.h"
#include "GeneDBDuplicates.h"
#include "GeneDBLocus.h"
#include <map>
#include <string>
#include <sstream>

using namespace std;

class GeneDBChromosome {
 public:
  string chromosomeName;

  // Data storage
  vector<GeneDBExon> exons;
  // Indices into the data
  
  map<int,int> byStart;
  multimap<string,int> byGene;
  multimap<string,int> byGeneDuplicates;

  bool LookupIndexByStart(int start, int &index) {
    map<int,int>::iterator it;
    it = byStart.find(start);
    if (it != byStart.end()) {
      index = (*it).second;
      return true;
    }
    else {
      return false;
    }
  }
  

  GeneDBDuplicates *duplicates;

  GeneDBChromosome() {
    duplicates = NULL;
  }

  void BuildIndex() {
    int i;
    for (i = 0; i < exons.size(); i++) {
      byStart[exons[i].start] = i;
      byGene.insert(pair<string,int>(exons[i].gene, i));
      byGeneDuplicates.insert(pair<string,int>(exons[i].geneDuplicates, i));
    }
  }

  void Read(istream &in) {
    string line;
    while(in and line != "START_CHROMOSOME") {
      getline(in,line);
      if (line == "START_CHROMOSOME") {
        continue;
      }
      stringstream strm(line);
      GeneDBExon exon;
      strm >> exon.gene >> exon.order >> exon.start >> exon.end >> exon.geneIndex >> exon.geneDuplicates;

      if (duplicates != NULL and 
          exon.geneDuplicates != "NONE" and
          exon.order == 0) {

        GeneDBLocus locus;
        locus.geneName      = exon.gene;
        locus.chromosome    = chromosomeName;
        locus.startPosition = exon.start;
        duplicates->AddDuplicate(exon.geneDuplicates, locus);
      }
      exons.push_back(exon);
    }
  }

  void Print(ostream &out) {



  }
  int size() {
    return exons.size();
  }

};


#endif
