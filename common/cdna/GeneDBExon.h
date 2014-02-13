#ifndef CDNA_GENE_DB_EXON_H_
#define CDNA_GENE_DB_EXON_H_

#include <string>
#include <iostream>

using namespace std;
class GeneDBExon {
 public:
  string gene;
  string geneDuplicates;
  int start, end;
  int order;
  string chromosome;
  int geneIndex;
  void Print(ostream &out) const {
    out << gene << " " << geneDuplicates << " " << start << " " << end << " " << order << " " << chromosome << endl;
  }

};

#endif
