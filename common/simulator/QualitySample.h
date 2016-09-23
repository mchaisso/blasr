#ifndef SIMULATOR_QUALITY_SAMPLE
#define SIMULATOR_QUALITY_SAMPLE
#include <vector>
#include <iostream>
#include "SMRTSequence.h"
#include "Types.h"
#include "qvs/QualityValue.h"
using namespace std;
#define NQV 5
#define NFV 3
#define NT  2
class QualitySample {
	public:
	QualityValue qv[NQV];
	HalfWord     frameValues[NFV];
	Nucleotide   tags[NT];

  void CopyFromSequence(SMRTSequence &seq, int pos) {
    qv[0] = seq.qual[pos];
    qv[1] = seq.deletionQV[pos];
    qv[2] = seq.insertionQV[pos];
    qv[3] = seq.substitutionQV[pos];
		qv[4] = seq.mergeQV[pos];
		assert(qv[3] != 0);
		
    tags[0] = seq.deletionTag[pos];
    tags[1] = seq.substitutionTag[pos];
    frameValues[0] = seq.pulseIndex[pos];
    frameValues[1] = seq.preBaseFrames[pos];
    frameValues[2] = seq.widthInFrames[pos];
  }

	void Write(ofstream &out) {
		out.write((char*)qv, sizeof(qv[0])*NQV);
		out.write((char*)frameValues, sizeof(frameValues[0])*NFV);
		out.write((char*)tags, sizeof(tags[0])*NT);
	}

	void Read(ifstream &in) {
		in.read((char*)qv, sizeof(QualityValue)*NQV);
		in.read((char*)frameValues, sizeof(frameValues[0])*NFV);
		in.read((char*)tags, sizeof(tags[0])*NT);
	}

	void CreateFromRead(SMRTSequence &seq, DNALength pos) {
		/*
		 * For now create with a static set of values:
		 * qvs:
		 * 0 QualityValue
		 * 1 DeletionQV
		 * 2 InsertionQV
		 * 3 SubstitutionQV
		 * 4 MergeQV
		 * tags:
		 * 0 DeletionTag
		 * 1 SubstitutionTag
		 * PulseData
		 * 0 PulseIndex
		 * 1 PreBaseFrames
		 * 2 WidthInFrames
		 */
		//		qv.resize(4);
		fill(&qv[0], &qv[NQV], 0);
		if (seq.qual.Empty() == false){ 
			qv[0] = seq.qual[pos];
		}
		if (seq.deletionQV.Empty() == false) {
			qv[1] = seq.deletionQV[pos];
		}
		if (seq.insertionQV.Empty() == false ){ 
			qv[2] = seq.insertionQV[pos];
		}
		if (seq.substitutionQV.Empty() == false){ 
			qv[3] = seq.substitutionQV[pos];
		}
		if (seq.mergeQV.Empty() == false){ 
			qv[4] = seq.mergeQV[pos];
		}
		fill(&tags[0], &tags[NT], 0);
		if (seq.deletionTag){ 
			tags[0] = seq.deletionTag[pos];
		}
		if (seq.substitutionTag) {
			tags[1] = seq.substitutionTag[pos];
		}
		
		std::fill(&frameValues[0], &frameValues[NFV], 0);
		if (seq.pulseIndex) {
			frameValues[0] = seq.pulseIndex[pos];
		}

		if (seq.preBaseFrames) {
			frameValues[1] = seq.preBaseFrames[pos];
		}

		if (seq.widthInFrames) {
			frameValues[2] = seq.widthInFrames[pos];
		}
	}
};


#endif
