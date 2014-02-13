#ifndef SIMULATOR_OUTPUT_SAMPLE_H_
#define SIMULATOR_OUTPUT_SAMPLE_H_

#include "QualitySample.h"
#include "SMRTSequence.h"

class OutputSample {
 public:
  enum Type {Match, Insertion, Deletion, Substitution};
  
  vector<QualitySample> qualities;
  vector<Nucleotide>    nucleotides;


  void Resize(int size) {
    qualities.resize(size);
    nucleotides.resize(size);
  }
  Type type;
  
  int CopyFromSeq(SMRTSequence &seq, int pos, int length=1) {
    Resize(length);
    int i;
    for (i = 0; i < length; i++) {
      qualities[i].CopyFromSequence(seq, pos+i);
      nucleotides[i] = seq.seq[pos+i];
    }
  }
  
  void Write(ofstream &out) {

    out.write((char*) &type, sizeof(type));
    int nNuc = nucleotides.size();

    out.write((char*)&nNuc, sizeof(int));
    int i;
    for (i = 0; i < qualities.size(); i++) {
      qualities[i].Write(out);
    }
    assert(nNuc == qualities.size());
    out.write((char*) &nucleotides[0], sizeof(Nucleotide)*nucleotides.size());

  }

  void Read(ifstream &in) {
    in.read((char*) &type, sizeof(Type));
    int nNuc;
    in.read((char*) &nNuc, sizeof(int));
    qualities.resize(nNuc);
    int i;
    for (i = 0; i < nNuc; i++) {
      qualities[i].Read(in);
    }
    nucleotides.resize(nNuc);
    in.read((char*) &nucleotides[0], sizeof(Nucleotide)* nNuc);

  }
};

#endif
