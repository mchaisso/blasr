#ifndef SIMULATOR_OUTPUT_SAMPLE_H_
#define SIMULATOR_OUTPUT_SAMPLE_H_

#include "QualitySample.h"
#include "SMRTSequence.h"

class OutputSample {
 public:
  enum Type {Match, Insertion, Deletion, Substitution, Merge};
  
  vector<QualitySample> qualities;
  vector<Nucleotide>    nucleotides;
	int nNuc;

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
			assert(nucleotides[i] >= 'A' and nucleotides[i] <= 'T');
    }
		nNuc = length;
  }
  
  void Write(ofstream &out) {
		
    out.write((char*) &type, sizeof(Type));
    out.write((char*) &nNuc, sizeof(int));
    int i;
		if (type != Deletion) {
			if (nNuc > 0) {
				for (i = 0; i < qualities.size(); i++) {
					qualities[i].Write(out);
				}
				out.write((char*) &nucleotides[0], sizeof(Nucleotide)*nucleotides.size());
			}
		}
  }

  void Read(ifstream &in) {
    in.read((char*) &type, sizeof(Type));
    in.read((char*) &nNuc, sizeof(int));
		if (type != Deletion) {
			qualities.resize(nNuc);
			int i;
			for (i = 0; i < nNuc; i++) {
				qualities[i].Read(in);
			}
			nucleotides.resize(nNuc);
			if (nNuc > 0) {
				in.read((char*) &nucleotides[0], sizeof(Nucleotide)* nNuc);
			}
		}
  }
};

#endif
