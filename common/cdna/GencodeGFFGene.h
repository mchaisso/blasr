#ifndef CDNA_GENCODE_GFF_GENE_H_
#define CDNA_GENCODE_GFF_GENE_H_

#include "GencodeGFFFile.h"
#include "FASTASequence.h"
#include "statistics/statutils.h"

class CompareByStartPosition {
 public:
  bool operator()(GencodeGFFEntry *a, GencodeGFFEntry *b) {
    return a->start < b->start;
  }
};

class ExonCoordinates {
 public:
  int start, end;
 ExonCoordinates(int s, int e) : start(s), end(e) {}

};

class GeneCoordinates {
 public:
  vector<ExonCoordinates> exonCoordinates;
  string chromosome;
  char strand;
};

class GencodeGFFGene {
 public:
  vector<GencodeGFFEntry*> exons;
  string geneName;
  string chromosome;
  char   strand;
  unsigned int GetStart();
  char GetStrand() {
    if (exons.size() == 0) {
      return 0;
    }
    else {
      return exons[0]->strand;
    }
  }

  void OrderExonsByStart() {
    CompareByStartPosition compare;
    sort(exons.begin(), exons.end(), compare);
  }

  void RecursiveGenerateSequence(FASTASequence &chromosome, 
                                 vector<FASTASequence> &sequences, 
                                 vector<GeneCoordinates> &sequenceCoordinates, 
                                 int curExon, string curSequence, GeneCoordinates &curGeneCoordinates,
                                 bool spliceRandomly=false) {
    if (curExon >= exons.size()) {
      //
      // Processed an entire sequence, add that to the list.
      //

      // Do nothing on empty sequences
      if (curSequence.size() == 0) {
        return;
      }
      
      FASTASequence geneSeq;
      
      geneSeq.Allocate(curSequence.size());
      memcpy((char*) geneSeq.seq, curSequence.c_str(), curSequence.size());
      geneSeq.length = curSequence.size();
      stringstream titleStrm;
      titleStrm << geneName;
      if (spliceRandomly) {
        titleStrm << "." << sequences.size();
      }
      geneSeq.CopyTitle(titleStrm.str());
      sequences.push_back(geneSeq);
      curGeneCoordinates.chromosome = this->chromosome;
      curGeneCoordinates.strand     = GetStrand();
      sequenceCoordinates.push_back(curGeneCoordinates);
    }
    else {
      //
      // Add the current exon to the sequence.
      //
      string sequenceWithCurExon = curSequence;
      int nextExon = curExon + 1;
      //
      // For sure, make a sequence with the following exon.
      //
      while (nextExon < exons.size() and 
             //
             // These two exons overlap
             //
             ((exons[nextExon]->start >= exons[curExon]->start and
               exons[nextExon]->start <= exons[curExon]->end) or
              (exons[nextExon]->end >= exons[curExon]->start and
               exons[nextExon]->end <= exons[curExon]->end) or 
              (exons[nextExon]->start <= exons[curExon]->start and
               exons[nextExon]->end >= exons[curExon]->end))) {
        nextExon = nextExon + 1;
      }
      int includedExon = curExon;
      if (spliceRandomly) {
        int randomExonOffset = RandomInt(nextExon - curExon);
        includedExon += randomExonOffset;
      }
      
      sequenceWithCurExon.append((char*)&chromosome.seq[exons[includedExon]->start],
                                 exons[includedExon]->end - exons[includedExon]->start);

      curGeneCoordinates.exonCoordinates.push_back(ExonCoordinates(exons[includedExon]->start, exons[includedExon]->end));

      //
      // Generate a sequence with the current exon
      //
      RecursiveGenerateSequence(chromosome, 
                                sequences, sequenceCoordinates,
                                nextExon, sequenceWithCurExon, curGeneCoordinates, spliceRandomly);

    }
  }

  void GenerateGeneSequences(FASTASequence &chromosome, vector<FASTASequence> &sequences, vector<GeneCoordinates> &coordinates,
                             bool generateAll=false) {
    string blank = "";
    GeneCoordinates emptyGeneCoordinates;
    RecursiveGenerateSequence(chromosome, sequences, coordinates, 0, blank, emptyGeneCoordinates, generateAll);
  }

};

void IndexGencodeGenes(GencodeGFFFile &gencodeFile, vector<GencodeGFFGene> &genes, string geneType = "protein_coding") {
  string curGeneName = "";
  int i;
  if (gencodeFile.entries.size() == 0) {
    return;
  }
  int curGene = 0;
  for (i = 0; i < gencodeFile.entries.size(); i++) {
    if (gencodeFile.entries[i].geneType == geneType) {
      if (curGeneName != gencodeFile.entries[i].geneName) {
        genes.push_back(GencodeGFFGene());
        curGene = genes.size() - 1;
        curGeneName = genes[curGene].geneName = gencodeFile.entries[i].geneName;
        genes[curGene].chromosome  = gencodeFile.entries[i].chr;
        genes[curGene].strand      = genes[curGene].GetStrand();
      }
      if (gencodeFile.entries[i].genomicLocusType == "exon") {
        genes[curGene].exons.push_back(&gencodeFile.entries[i]);
      }
    }
  }        
}
  

#endif
