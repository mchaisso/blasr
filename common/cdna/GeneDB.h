#ifndef CDNA_GENE_DB_H_
#define CDNA_GENE_DB_H_

#include "GeneDBChromosome.h"
#include "GeneDBDuplicates.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <map>

class ExonIndex {
 public:
  string chromosome;
  int    index;
  ExonIndex(string c, int i) {
    chromosome = c;
    index = i;
  }
};

typedef map<string, GeneDBChromosome*> MapGeneDBChromosome;
using namespace std;
class GeneDB : public MapGeneDBChromosome {
 public:
  multimap<string, ExonIndex> geneToExon;

  vector<GeneDBChromosome*> chromosomes;
  GeneDBDuplicates geneDBDuplicates;
  GeneDBChromosome* Find(string key) {
    MapGeneDBChromosome::iterator it = find(key);
    if ( it == end() ) {
      return NULL;
    }
    else {
      return it->second;
    }
  }

  GeneDBChromosome& operator[](string key) {
    MapGeneDBChromosome::iterator it = find(key);
    if ( it == end() ) {
      GeneDBChromosome *newChromosome = AddChromosome(key);
      return *newChromosome;
    }
    else {
      return *it->second;
    }
  }
  
  GeneDBChromosome* AddChromosome(string chrName) {

    MapGeneDBChromosome::iterator it = find(chrName);
    if ( it != end() ) {
      return it->second;
    }
    GeneDBChromosome* newChromosome = new GeneDBChromosome;
    insert(pair<string,GeneDBChromosome*>(chrName, newChromosome));
    newChromosome->chromosomeName = chrName;

    // 
    // Link to the duplicates structure so these may be stored as the chromosome is read in.
    //
    newChromosome->duplicates = &geneDBDuplicates;
    return newChromosome;
  }

  void IndexChromosomes() {
    GeneDB::iterator it = begin(), endIt = end();
    for (; it != endIt; ++it) {
      (*it).second->BuildIndex();
      int e;
      for (e = 0; e < (*it).second->exons.size(); e++) {
        if ((*it).second->exons[e].geneDuplicates != "NONE") {
          geneToExon.insert(pair<string, ExonIndex>((*it).second->exons[e].gene, ExonIndex((*it).first, e)));
        }
      }
    }
  }
    
 bool Read(string inFileName) {
    ifstream in;
    CrucialOpen(inFileName, in, std::ios::in);
    return Read(in);
  }

  bool Read(istream &in) {
    string line;
    getline(in, line);
    if (line != "START_CHROMOSOME") {
      cout << "Malformatted GeneDB file." << endl;
      return false;
    }
    while(in) {
      string chrName;
      in >> chrName;
      GeneDBChromosome *chromosome;
      chromosome = AddChromosome(chrName);
      chromosome->Read(in);
    }
  }

  void Write(string outFileName) {
    ofstream outFile;
    CrucialOpen(outFileName, outFile, std::ios::out);
    Write(outFile);
  }

  void Write(ostream &out) {
    MapGeneDBChromosome::iterator it = begin();    
    MapGeneDBChromosome::iterator endIt = end();
    while(it != endIt) {
      if ((*it).second->exons.size() == 0) {
        ++it;
        continue;
      }
      out << "START_CHROMOSOME" << endl;
      out << (*it).second->chromosomeName << endl;
      (*it).second->Print(cout);
    }
  }
};


#endif
