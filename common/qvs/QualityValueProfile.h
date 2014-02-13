#ifndef QVS_QUALITY_VALUE_PROFILE_H_
#define QVS_QUALITY_VALUE_PROFILE_H_

#include "../tuples/DNATuple.h"
#include "../tuples/TupleMetrics.h"
#include "../datastructures/matrix/FlatMatrix.h"
#include "QualityValue.h"

class QualityValueProfile {
	int wordSize;
	int numQualityValues;
	FlatMatrix2D<int> profile;
	int nWords;
	TupleMetrics tm;
 public:
	static const int CDF_GRANULARITY = 10000;
	QualityValueProfile(int _wordSize, int _numQualityValues) {
		wordSize = _wordSize;
		numQualityValues = _numQualityValues;
		// Initialize the tuple metrics to map from seq->index
		tm.Initialize(wordSize);
		nWords = 1 << (2*wordSize);
		// Initialize the matrix of quality values.
		profile.Grow(nWords, numQualityValues);
		profile.Initialize(0);
	}

	void Update(Nucleotide *seq, QualityValue qv) {
		DNATuple tuple;
		if (tuple.FromStringLR(seq, tm)) {
			profile.Set(tuple.tuple, qv, profile(tuple.tuple, qv) + 1);
		}
	}

	void Print(ofstream &out) {
		out << wordSize << " " << numQualityValues << " " << CDF_GRANULARITY << endl;
		profile.Print(out);
	}

	void ProfileToCDF() {
		int qv;
		int wordIndex;
		for (wordIndex = 0 ; wordIndex < nWords; wordIndex++) {
			int totalSamples = 0;			
			for (qv = 0; qv < numQualityValues; qv++) {
				totalSamples += profile(wordIndex, qv);
				profile.Set(wordIndex, qv, totalSamples);
			}
			for (qv = 0; qv < numQualityValues; qv++ ){ 
				profile.Set(wordIndex, qv, ((profile(wordIndex,qv) * 1.0) / totalSamples)  * CDF_GRANULARITY);
			}
		}
	}
};


#endif
