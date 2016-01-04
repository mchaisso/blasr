#ifndef DATA_HDF_HDF_SCAN_DATA_WRITER_H_
#define DATA_HDF_HDF_SCAN_DATA_WRITER_H_

#include "HDFGroup.h"
#include "HDFAtom.h"
#include "PlatformId.h"
#include "Enumerations.h"
#include "datastructures/reads/ScanData.h"
//
// The SanDataWriter cannot live outside 

class HDFScanDataWriter {
 public:
	bool fileHasScanData, useRunCode, fileHasBaseMap, fileHasDyeSet;
	HDFGroup scanDataGroup;
	HDFGroup dyeSetGroup;
	HDFGroup acqParamsGroup;
	HDFGroup runInfoGroup;
	bool initializedAcqParamsGroup, initializedRunInfoGroup;
	bool useWhenStarted;
	HDFAtom<string> whenStartedAtom;
	HDFAtom<unsigned int> platformIdAtom;
	HDFAtom<float> frameRateAtom;
	HDFAtom<unsigned int> numFramesAtom;
	HDFAtom<string> movieNameAtom;
	HDFAtom<string> runCodeAtom;
	HDFAtom<string> baseMapAtom;
	HDFAtom<string> sequencingKitAtom;
	HDFAtom<string> sequencingChemistryAtom;
	HDFAtom<string> basecallerVersionAtom;

	//
	// It is useful to cache the movie name in the reader since this is
	// loaded once upon initialization, and may be fetched when loading
	// reads one at a time.
	//
	bool   useMovieName;
	string movieName, runCode;
	map<char, int> baseMap;
	PlatformId platformId;

	HDFScanDataWriter() {
		//
		// Assume the file is written without a movie name.  This is
		// flipped when a movie name is found.
		//
		useMovieName    = false;
		useRunCode      = false;
		useWhenStarted  = false;
		fileHasScanData = false;
		movieName = "";
		runCode   = "";
		platformId      = NoPlatform;
		initializedAcqParamsGroup = initializedRunInfoGroup = false;
	}

	void Initialize(HDFGroup &parentGroup) {

		parentGroup.AddGroup("ScanData"); 

		if (scanDataGroup.Initialize(parentGroup, "ScanData") == 0) {
      cout << "ERROR, could not create ScanData." << endl;
      exit(1);
    }

    scanDataGroup.AddGroup("RunInfo");
    if (runInfoGroup.Initialize(scanDataGroup, "RunInfo") == 0) {
      cout << "ERROR, could not create creating /ScanData/RunInfo"<< endl;
      exit(1);
    }

    scanDataGroup.AddGroup("AcqParams");
    if (acqParamsGroup.Initialize(scanDataGroup, "AcqParams") == 0) {
      cout << "ERROR creating hdf acq params." << endl;
      exit(1);
    }
	}
	
	void Write(ScanData &scanData) {

		/*
			items to write
 public:
	PlatformId platformId;
	float frameRate;
	unsigned int numFrames;
	string movieName, runCode;
	string whenStarted;
    map<char, int> baseMap;
	string GetMovieName() {
		return movieName;
	}
		*/

		


	}

};

#endif
