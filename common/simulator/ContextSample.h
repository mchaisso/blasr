#ifndef SIMULATOR_CONTEXT_SAMPLE_H_
#define SIMULATOR_CONTEXT_SAMPLE_H_

#include<string>
#include<vector>
#include<iostream>

#include "simulator/QualitySample.h"
#include "statistics/statutils.h"
class ContextSample {
 public:
	vector<QualitySample> samples;
	int minSamples;
	int maxSamples;
	int reachedMinSamples;
	ContextSample() {
		minSamples = 0;
		maxSamples = 0;
		reachedMinSamples = 0;
	}
	int GetNSamples() {
		return samples.size();
	}
	int AppendSample(SMRTSequence &seq, DNALength pos) {
		// 
		// Check to see if the space has been sampled enough.
		//
		if (maxSamples != 0 and samples.size() >= maxSamples)
			return 0;

		samples.resize(samples.size()+1);
		samples[samples.size()-1].CreateFromRead(seq,pos);
		if (samples.size() > minSamples and reachedMinSamples == 0) {
			reachedMinSamples = 1;
			return 1;
		}
		else {
			return 0;
		}
	}

	QualitySample* GetRandomQualitySample() {
		int sampleIndex = RandomInt(samples.size());
		return &samples[sampleIndex];
	}

	void Write(ofstream &out) {
		int s;
		int sampleSize = samples.size();
		out.write((char*) &sampleSize, sizeof(sampleSize));
		for (s = 0; s < samples.size(); s++ ){
			samples[s].Write(out);
		}
	}
	void Read(ifstream &in) {
		int s;
		int sampleSize;
		in.read((char*)&sampleSize, sizeof(sampleSize));
		samples.resize(sampleSize);
		for (s = 0; s < sampleSize; s++ ){
			samples[s].Read(in);
		}
	}
};

#endif

