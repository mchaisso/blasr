#ifndef DATASTRUCTURES_READS_SCAN_DATA_H_
#define DATASTRUCTURES_READS_SCAN_DATA_H_

#include "../../Enumerations.h"
#include "../../data/hdf/PlatformId.h"
class ScanData {
 public:
	PlatformId platformId;
	float frameRate;
	unsigned int numFrames;
	string movieName, runCode;
	string bindingKit;
	string whenStarted;
	string sequencingChemistry;
	string sequencingKit;
	string basecallerVersion;
    map<char, int> baseMap;
	string GetMovieName() {
		return movieName;
	}
	void Clear() {
		frameRate = 0;
		numFrames = 0;
		movieName = 
			runCode = 
			bindingKit = 
			sequencingChemistry =
			sequencingKit =
			basecallerVersion = 
			whenStarted = "";

	}

};

#endif
