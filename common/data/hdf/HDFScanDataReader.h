#ifndef DATA_HDF_HDF_SCAN_DATA_READER_H_
#define DATA_HDF_HDF_SCAN_DATA_READER_H_

#include "HDFGroup.h"
#include "HDFAtom.h"
#include "PlatformId.h"
#include "Enumerations.h"
#include "datastructures/reads/ScanData.h"
//
// The SanDataReader cannot live outside 

class HDFScanDataReader {
 public:
	bool fileHasScanData, 
		useRunCode, 
		fileHasBaseMap,
		fileHasDyeSet,
		useBindingKit,
		useSequencingKit,
		useSequencingChemistry,
		useBasecallerVersion;
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
	HDFAtom<string> bindingKitAtom;
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
	HDFScanDataReader() {
		//
		// Assume the file is written without a movie name.  This is
		// flipped when a movie name is found.
		//
		useMovieName    = false;
		useRunCode      = false;
		useWhenStarted  = false;
		fileHasScanData = false;
		useBindingKit   = false;
		useSequencingKit = false;
		useSequencingChemistry=false;
		useBasecallerVersion=false;

		movieName = "";
		runCode   = "";
		platformId      = NoPlatform;
		initializedAcqParamsGroup = initializedRunInfoGroup = false;
	}

	int InitializeAcqParamsAtoms() {
		if (frameRateAtom.Initialize(acqParamsGroup.group, "FrameRate") == 0) { return 0; }
		if (numFramesAtom.Initialize(acqParamsGroup.group, "NumFrames") == 0) { return 0; }
		if (acqParamsGroup.ContainsAttribute("WhenStarted")) {
			if (whenStartedAtom.Initialize(acqParamsGroup.group, "WhenStarted") == 0) { return 0; }
			useWhenStarted = true;
		}
		return 1;
	}

	//
	// This is created on top of a file that is already opened, so
	// instead of initializing by opening a file, it is initialized by
	// passing the root group that should contain the ScanData group.  
	// When shutting down, no file handle needs to be closed.
	//
	int Initialize(HDFGroup *pulseDataGroup) {

		//
		// Initiailze groups for reading data.
		//
		
		initializedAcqParamsGroup = false;
		initializedRunInfoGroup   = false;
		if (pulseDataGroup->ContainsObject("ScanData") == 0 or 
				scanDataGroup.Initialize(pulseDataGroup->group, "ScanData") == 0) {
			//			return 0;
			fileHasScanData = false;
		}
		else {
			fileHasScanData = true;
		}

		if (scanDataGroup.ContainsObject("DyeSet") == 0 or
				dyeSetGroup.Initialize(scanDataGroup.group, "DyeSet") == 0) {
			fileHasDyeSet = false;
			//			return 0;
		}
		else {
			fileHasDyeSet = true;
		}
		if (scanDataGroup.ContainsObject("AcqParams") == 0 or
		    acqParamsGroup.Initialize(scanDataGroup.group, "AcqParams") == 0) {
			initializedAcqParamsGroup = false;
			//			return 0;
		}
		else {
			initializedAcqParamsGroup = true;
		}

		if (scanDataGroup.ContainsObject("RunInfo") == 0 or
				runInfoGroup.Initialize(scanDataGroup.group, "RunInfo") == 0) {
			//			return 0;
			initializedRunInfoGroup = false;
		}
		else {
			initializedRunInfoGroup = true;
		}
		if (InitializeAcqParamsAtoms() == 0) {
			//			return 0;
		}

		//
		// Read in the data that will be used later on either per read or
		// when the entire bas/pls file is read.
		//
		
		if (ReadPlatformId(platformId) == 0) {
			//			return 0;
		}

		if (runInfoGroup.ContainsAttribute("RunCode") and
				runCodeAtom.Initialize(runInfoGroup, "RunCode")) {
			useRunCode = true;
		}

		if (runInfoGroup.ContainsAttribute("BindingKit") and
				bindingKitAtom.Initialize(runInfoGroup, "BindingKit")) {
			useBindingKit = true;
		}

		if (runInfoGroup.ContainsAttribute("SequencingKit") and
				sequencingKitAtom.Initialize(runInfoGroup, "SequencingKit")) {
			useSequencingKit = true;
		}

		if (runInfoGroup.ContainsAttribute("SequencingChemistry") and
				sequencingChemistryAtom.Initialize(runInfoGroup, "SequencingChemistry")) {
			useSequencingChemistry = true;
		}


		//
		// Load baseMap which maps bases (ATGC) to channel orders.
		// This should always be present.
		//
		if (fileHasDyeSet) {
			if (LoadBaseMap(baseMap) == 0) {
				return 0;
			}
		}

		//
		// Attempt to load the movie name.  This is not always present.
		//
		LoadMovieName(movieName);

		return 1;
	}
	
	string GetRunCode() {
		if (useRunCode) {
			string runCode;
			runCodeAtom.Read(runCode);
			return runCode;
		}
		else {
			return "";
		}
	}
	
	string GetMovieName() {
		LoadMovieName(movieName);
		return movieName;
	}
	
	int Read(ScanData &scanData) {
		// All parameters below are required.
		if (ReadPlatformId(scanData.platformId) == 0) return 0;
		LoadMovieName(scanData.movieName);
		LoadBaseMap(scanData.baseMap);
		
		if (useRunCode) {
			runCodeAtom.Read(scanData.runCode);
		}
		if (useBindingKit) {
			bindingKitAtom.Read(scanData.bindingKit);
		}
		if (useSequencingChemistry) {
			bindingKitAtom.Read(scanData.sequencingChemistry);
		}

		if (useSequencingKit) {
			sequencingKitAtom.Read(scanData.sequencingKit);
		}

		if (useBasecallerVersion) {
			basecallerVersionAtom.Read(scanData.basecallerVersion);
		}

		frameRateAtom.Read(scanData.frameRate);
		numFramesAtom.Read(scanData.numFrames);
		
		if (useWhenStarted) {
			whenStartedAtom.Read(scanData.whenStarted);
		}

	}

	void ReadWhenStarted(string &whenStarted) {
		whenStartedAtom.Read(whenStarted);
	}

	PlatformId GetPlatformId() {
		return platformId;
	}

	int ReadPlatformId(PlatformId &pid) {
		if (runInfoGroup.ContainsAttribute("PlatformId")) {
			if (platformIdAtom.Initialize(runInfoGroup, "PlatformId") == 0) {
				return 0;
			}
			platformIdAtom.Read((unsigned int&)pid);
		}
		else {
			pid = AstroPlatform;
		}
		return 1;
	}

	int LoadMovieName(string &movieName) {
		// Groups for building read names
		if (runInfoGroup.ContainsAttribute("MovieName") and
				movieNameAtom.Initialize(runInfoGroup, "MovieName")) {
			useMovieName = true;
			movieNameAtom.Read(movieName);
			int e = movieName.size() - 1;
			while (e > 0 and movieName[e] == ' ') e--;
			movieName= movieName.substr(0,e+1);
			return 1;
		}
		else {
			return 0;
		}
	}

	int LoadBaseMap(map<char, int> & baseMap) {
		// Map bases to channel order in hdf pls file.
		if (fileHasDyeSet) {
			if (dyeSetGroup.ContainsAttribute("BaseMap") and
					baseMapAtom.Initialize(dyeSetGroup, "BaseMap")) {
				string baseMapStr;
				baseMapAtom.Read(baseMapStr);
				if (baseMapStr.size() != 4) {
					cout << "ERROR, there are more than four types of bases "
							 << "according to /ScanData/DyeSet/BaseMap." << endl;
					exit(1);
				}
				baseMap.clear();
				for(int i = 0; i < baseMapStr.size(); i++) {
					baseMap[toupper(baseMapStr[i])] = i;
					baseMap[tolower(baseMapStr[i])] = i;
				}
				return 1;
			}
		}
		return 0;
	}

	void Close() {
		if (useMovieName) {
			movieNameAtom.dataspace.close();
		}
		if (useRunCode) {
			runCodeAtom.dataspace.close();
		}
		if (useWhenStarted) {
			whenStartedAtom.dataspace.close();
		}
		baseMapAtom.dataspace.close();
		platformIdAtom.dataspace.close();
		frameRateAtom.dataspace.close();
		numFramesAtom.dataspace.close();

		scanDataGroup.Close();
		dyeSetGroup.Close();
		acqParamsGroup.Close();
		runInfoGroup.Close();
	}


};

#endif
