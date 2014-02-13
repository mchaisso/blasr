#include "../common/files/ReaderAgglomerate.h"
#include "../common/SMRTSequence.h"
#include "../common/utils/FileOfFileNames.h"
#include "../common/simulator/ContextSet.h"

void PrintUsage() {
	cout << "storeQualityByContext - grab quality values from .bas.h5 files until minimum requirements for the number of times a context has been sampled are met." << endl;
	cout << "usage: storeQualityByContext bas.h5|fofn  output.qbc  [options] " << endl;
	cout << "options: " << endl
			 << " -contextlength L  The length of the context to sample" << endl
			 << " -minSamples S(500)Report pass if all contexts are sampled" <<endl
			 << "                   at least S times." << endl
			 << " -maxSamples S(1000)Stop sampling a context once it has reached" << endl
			 << "                   S samples." << endl
			 << " -minAvgQual q     Only sample a read if it has a minimum average quality value of q." << endl;
}

int main(int argc, char* argv[]) {
	string inFileName, outFileName;
	int contextLength = 5;
	int minSamples = 500;
	int maxSamples = 1000;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	
	int argi = 1;
	inFileName=  argv[argi++];
	outFileName = argv[argi++];
	int minAverageQual = 0;

	while (argi < argc) {
		if (strcmp(argv[argi], "-contextLength") == 0) {
			contextLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minSamples") == 0) {
			minSamples = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-maxSamples") == 0) {
			maxSamples = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minAvgQual") == 0) {
			minAverageQual = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			cout << "Bad option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}

	SMRTSequence read;
	vector<string> inputFileNames;
	FileOfFileNames::StoreFileOrFileList(inFileName, inputFileNames);
	ofstream sampleOut;
	CrucialOpen(outFileName, sampleOut, std::ios::out|std::ios::binary);
	int fileNameIndex;
	ContextSampleMap samples;
	
	int numContextsReached = 0;
	int numContexts = 1 << (contextLength*2);
	ReaderAgglomerate reader;
	samples.contextLength = contextLength;
	for (fileNameIndex = 0; 
			 fileNameIndex < inputFileNames.size() and
				 numContextsReached < numContexts ; fileNameIndex++) {
		reader.Initialize(inputFileNames[fileNameIndex]);
		while(reader.GetNext(read) and numContextsReached < numContexts) {
			if (read.length < contextLength) {
				continue;
			}
			int totalQual = 0;
			int readPos;
			for (readPos = 0; readPos < read.length - contextLength + 1 ; readPos++) {			
				totalQual += read.qual[readPos];
			}
			if (totalQual / (read.length - contextLength + 1) < minAverageQual)
				continue;
			for (readPos = 0; readPos < read.length - contextLength + 1 and numContextsReached < numContexts; readPos++) {
				string context;
				context.assign((const char*) &read.seq[readPos],contextLength);
				// Make sure this context exists.
				if (samples.find(context) == samples.end()) {
					samples[context] = new ContextSample;
					samples[context]->minSamples = minSamples;
					samples[context]->maxSamples = maxSamples;
				}
				if (samples[context]->AppendSample(read, readPos) == 1) {
					++numContextsReached;
				}
			}
		}
		reader.Close();
	}
	cout  << numContextsReached << " of "  << numContexts << " sampled to " << minSamples << endl;
	samples.Write(sampleOut);

	return 0;
}

	
