#ifndef DATA_HDF_HDF_BAS_WRITER_H_
#define DATA_HDF_HDF_BAS_WRITER_H_

#include "HDFArray.h"
#include "BufferedHDFArray.h"
#include "HDF2DArray.h"
#include "BufferedHDF2DArray.h"
#include "HDFAtom.h"
#include "HDFFile.h"
#include "DatasetCollection.h"
#include "Enumerations.h"
#include "utils/SMRTReadUtils.h"
#include "FASTQSequence.h"
#include "Types.h"
#include <sstream>
#include "HDFScanDataWriter.h"

using namespace H5;
using namespace std;

class HDFBasWriter : public DatasetCollection {
	HDFFile outFile;
	string hdfFileName;
	string movieName, runCode;
	static const int bufferSize = 16;
	
	float frameRate;
	float numFrames;
	string changeListID;

	HDFScanDataWriter scanDataWriter;

	HDFAtom<string> changeListIDAtom;

	BufferedHDFArray<int> nElemArray;
	BufferedHDFArray<int> zmwXCoordArray;
	BufferedHDFArray<int> zmwYCoordArray;
	BufferedHDFArray<unsigned char> baseArray;
	BufferedHDFArray<unsigned char> qualArray;
	BufferedHDFArray<unsigned int> simulatedCoordinateArray;
	BufferedHDFArray<unsigned int> simulatedSequenceIndexArray;

	HDFAtom<string> movieNameAtom, runCodeAtom, platformNameAtom;
	HDFAtom<string> sequencingKitAtom, bindingKitAtom;
	HDFAtom<unsigned int> platformIdAtom;

	//
	// Astro specific arrays.
	//
	BufferedHDF2DArray<int16_t> holeXY2D;

	//
	// Springfield specific arrays.
	//
	BufferedHDFArray<unsigned int> holeNumberArray;
	BufferedHDFArray<unsigned char> holeStatusArray;
	
	//
	// Define arrays for rich quality values.
	// 

	BufferedHDFArray<unsigned char> deletionQVArray;
	BufferedHDFArray<unsigned char> deletionTagArray;
	BufferedHDFArray<unsigned char> insertionQVArray;
	BufferedHDFArray<unsigned char> substitutionTagArray;
	BufferedHDFArray<unsigned char> substitutionQVArray;
	BufferedHDFArray<unsigned char> mergeQVArray;
	BufferedHDFArray<HalfWord> preBaseFramesArray;
	BufferedHDFArray<HalfWord> widthInFramesArray;
	BufferedHDFArray<int> pulseIndexArray;
	BufferedHDF2DArray<unsigned char> preBaseDeletionQVArray;
	HDFGroup rootGroup;

	HDFGroup baseCallGroup;
	HDFGroup zmwGroup;
	HDFGroup plsZMWGroup;
	PlatformType platformId;
	string   platformName;
 public:
	HDFGroup pulseDataGroup;
	~HDFBasWriter() { 
		//
		// Assume that flushing out and closing the hdf file must be one
		// manually and not in a destructor.
		//
	}
	void InitializeDefaultIncludedFields() {
		IncludeField("Basecall");
		IncludeField("DeletionQV");
		IncludeField("DeletionTag");
		IncludeField("InsertionQV");
		IncludeField("SubstitutionTag");
		IncludeField("SubstitutionQV");
		IncludeField("QualityValue");
		IncludeField("WidthInFrames");
		IncludeField("PulseIndex");
		IncludeField("PreBaseFrames");
		IncludeField("HoleNumber");    
    IncludeField("HoleStatus");
    IncludeField("MergeQV");
	}
	

	void Flush() {
		nElemArray.Flush();
		if (includedFields["zmwXCoord"])
		zmwXCoordArray.Flush();
		if (includedFields["zmwYCoord"])
			zmwYCoordArray.Flush();
		if (includedFields["Basecall"])
			baseArray.Flush();
		if (includedFields["QualityValue"])
			qualArray.Flush();
		if (includedFields["DeletionQV"])
			deletionQVArray.Flush();
		if (includedFields["DeletionTag"])
			deletionTagArray.Flush();
		if (includedFields["InsertionQV"])
			insertionQVArray.Flush();
		if (includedFields["SubstitutionTag"])
			substitutionTagArray.Flush();
		if (includedFields["SubstitutionQV"])
			substitutionQVArray.Flush();
		if (includedFields["HoleNumber"])
			holeNumberArray.Flush();
		if (includedFields["HoleStatus"])
			holeStatusArray.Flush();
		if (includedFields["PreBaseFrames"])
			preBaseFramesArray.Flush();
    if (includedFields["PulseIndex"]) 
      pulseIndexArray.Flush();
		if (includedFields["WidthInFrames"])
			widthInFramesArray.Flush();
		if (includedFields["HoleXY"])
			holeXY2D.Flush();
    if (includedFields["MergeQV"]) 
      mergeQVArray.Flush();
		if (includedFields["SimulatedCoordinate"])
			simulatedCoordinateArray.Flush();
		if (includedFields["SimulatedSequenceIndex"]) 
			simulatedSequenceIndexArray.Flush();
	}

	HDFBasWriter() {
		/*
		 * Default to astro for now.  This may need to change to a NO_ID
		 * platform, in which case it must be set with Initialize().
		 */
		frameRate = 100.0;
		numFrames = 10000000;
		fieldNames.push_back("zmwXCoord");
		fieldNames.push_back("zmwYCoord");
		fieldNames.push_back("QualityValue");
		fieldNames.push_back("Basecall");
		fieldNames.push_back("DeletionQV");
		fieldNames.push_back("DeletionTag");
		fieldNames.push_back("InsertionQV");
		fieldNames.push_back("SubstitutionQV");
		fieldNames.push_back("SubstitutionTag");
		fieldNames.push_back("MergeQV");
		fieldNames.push_back("WidthInFrames");
		fieldNames.push_back("HoleNumber");
		fieldNames.push_back("HoleStatus");
		fieldNames.push_back("HoleXY");
		fieldNames.push_back("PreBaseFrames");
		fieldNames.push_back("PulseIndex");
		fieldNames.push_back("SimulatedCoordinate");
		fieldNames.push_back("SimulatedSequenceIndex");
		InitializeAllFields(false);
		platformId = Springfield;
	}

	void Close() {
		Flush();
		outFile.Close();
	}

	void SetPlatform(PlatformType _platform) {
		platformId = _platform;
	}

	void SetChangeListID(string _changeListID) {
		changeListID = _changeListID;
	}

	void AddMovieName(string movieName) {
		movieNameAtom.Create(scanDataWriter.runInfoGroup.group, "MovieName",movieName);
	}
	
	/*
	 * Initialization without a runCode is implicitly a springfield
	 * platform.  You can change it if you really want.
	 */

	void Initialize(string _hdfFileName, string movieName, string _changeListID) {
		SetChangeListID(_changeListID);
		Initialize(_hdfFileName, Springfield);
		AddMovieName(movieName);
    AddPlatformInformation(Springfield);
	}

	void Initialize(string _hdfFileName, string movieName, PlatformType _platform = Springfield) {
		Initialize(_hdfFileName, _platform);
		AddMovieName(movieName);
    AddPlatformInformation(_platform);
	}

	void Initialize(string _hdfFileName, string movieName, string runCode, string _changeListID ) {
		Initialize(_hdfFileName, Springfield);
		if (movieName != "" and runCode != "") {
			AddRunInfo(movieName, runCode);
    }
    AddChangeListID(_changeListID);
    AddPlatformInformation(Springfield);
	}

  void AddPlatformInformation(PlatformType _platform) {
    platformId = _platform;
    AddPlatformId(platformId);
    AddPlatformName(platformId);
  }

	void AddChangeListID(string cl) {
    changeListIDAtom.Create(baseCallGroup.group, "ChangeListID", cl);
  }

	void AddRunInfo(string movieName, string runCode) {
		AddMovieName(movieName);
		runCodeAtom.Create(scanDataWriter.runInfoGroup.group, "RunCode", runCode);
	}

	void AddSequencingKit(string sequencingKit) {
		sequencingKitAtom.Create(scanDataWriter.runInfoGroup.group, "SequencingKit", sequencingKit);
	}

	void AddBindingKit(string bindingKit) {
		bindingKitAtom.Create(scanDataWriter.runInfoGroup.group, "BindingKit", bindingKit);
	}

	void AddPlatformName(PlatformType platformId) {
    if (platformId == Springfield) {
      platformNameAtom.Create(scanDataWriter.runInfoGroup.group, "PlatformName", "Springfield");
    }
    else if (platformId == Astro) {
      platformNameAtom.Create(scanDataWriter.runInfoGroup.group, "PlatformName", "Astro");
    }
	}

	void AddPlatformId(PlatformType _platformId) {
		platformIdAtom.Create(scanDataWriter.runInfoGroup.group, "PlatformId");
		platformIdAtom.Write((unsigned int)_platformId);
	}

	void WriteSimulatedCoordinate(unsigned int coord) {
		simulatedCoordinateArray.Write(&coord,1);
	}

	void WriteSimulatedSequenceIndex(unsigned int index) {
		simulatedSequenceIndexArray.Write(&index,1);
	}

	void Initialize(string _hdfFileName, PlatformType _platform,
            const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT) {
		hdfFileName = _hdfFileName;
		platformId  = _platform;
		//outFile.Open(hdfFileName, H5F_ACC_TRUNC);
		outFile.Open(hdfFileName, H5F_ACC_TRUNC, fileAccPropList);
		outFile.rootGroup.AddGroup("PulseData"); 

    if (pulseDataGroup.Initialize(outFile.rootGroup, "PulseData") == 0) {
      cout << "ERROR, could not create file " << _hdfFileName << ".  Error creating group /PulseData." << endl;
      exit(1);
    }
      
    pulseDataGroup.AddGroup("BaseCalls"); 
    if (baseCallGroup.Initialize(pulseDataGroup, "BaseCalls") == 0) {
      cout << "ERROR, could not create file " << _hdfFileName << ".  Error creating group /PulseData/BaseCall." << endl;
      exit(1);
    }

    
    baseCallGroup.AddGroup("ZMW"); 
    if (zmwGroup.Initialize(baseCallGroup, "ZMW") == 0) {
      cout << "ERROR, could not create file " << _hdfFileName << ".  Error creating group /PulseData/BaseCall/ZMW." << endl;
      exit(1);
    }

    
		scanDataWriter.Initialize(outFile.rootGroup);
      

		HDFAtom<float> frameRateAtom;
		HDFAtom<unsigned int> numFramesAtom;
    frameRateAtom.Create(scanDataWriter.acqParamsGroup.group, "FrameRate");
    numFramesAtom.Create(scanDataWriter.acqParamsGroup.group, "NumFrames");
    frameRateAtom.Write(frameRate);
    numFramesAtom.Write(numFrames);
		

		if (changeListID != "") {
			changeListIDAtom.Create(baseCallGroup.group, "ChangeListID", changeListID);
		}

		nElemArray.Initialize(zmwGroup, "NumEvent");
		if (includedFields["Basecall"])
			baseArray.Initialize(baseCallGroup, "Basecall");
		if (includedFields["QualityValue"])
			qualArray.Initialize(baseCallGroup, "QualityValue");
		if (includedFields["DeletionQV"])
			deletionQVArray.Initialize(baseCallGroup, "DeletionQV");
		if (includedFields["DeletionTag"])
			deletionTagArray.Initialize(baseCallGroup, "DeletionTag");
		if (includedFields["InsertionQV"])
			insertionQVArray.Initialize(baseCallGroup, "InsertionQV");
		if (includedFields["MergeQV"])
			mergeQVArray.Initialize(baseCallGroup, "MergeQV");
		if (includedFields["PreBaseDeletionQV"])
			preBaseDeletionQVArray.Initialize(baseCallGroup, "PreBaseDeletionQV", 4);
		if (includedFields["SubstitutionTag"])
			substitutionTagArray.Initialize(baseCallGroup,   "SubstitutionTag");
		if (includedFields["SubstitutionQV"])
			substitutionQVArray.Initialize(baseCallGroup,    "SubstitutionQV");
		if (includedFields["WidthInFrames"])
			widthInFramesArray.Initialize(baseCallGroup, "WidthInFrames");
		if (includedFields["PreBaseFrames"])
			preBaseFramesArray.Initialize(baseCallGroup, "PreBaseFrames");
		if (includedFields["SimulatedCoordinate"]) {
			simulatedCoordinateArray.Initialize(baseCallGroup, "SimulatedCoordinate");
		}
		if (includedFields["SimulatedSequenceIndex"]) {
			simulatedSequenceIndexArray.Initialize(baseCallGroup, "SimulatedSequenceIndex");
		}
    if (includedFields["PulseIndex"]) {
      pulseIndexArray.Initialize(baseCallGroup, "PulseIndex");
    }

		if (platformId == Astro or includedFields["HoleXY"]) {
			holeXY2D.Initialize(zmwGroup, "HoleXY", 2);
		}
    includedFields["HoleNumber"] = true;
		holeNumberArray.Initialize(zmwGroup, "HoleNumber");
    includedFields["HoleStatus"] = true;
		holeStatusArray.Initialize(zmwGroup, "HoleStatus");
	}

	int WriteHoleXY(int x=0, int y=0) {
		int16_t xy[2] = {(uint16_t) x, (uint16_t) y};
		holeXY2D.WriteRow(xy, 2);
	}		

	int WriteIdentifiers(UInt holeNumber, unsigned char holeStatus, int x=0, int y=0 ) {
		//
		// Write hole number regardless of platform type.
		//
		holeNumberArray.Write(&holeNumber, 1);
    holeStatusArray.Write(&holeStatus, 1);
		if (platformId == Astro or includedFields["HoleXY"]) {
			WriteHoleXY(x,y);
		}
		return 1;
	}
	

	int WriteQualities(FASTQSequence &seq) {
		qualArray.Write(seq.qual.data, seq.length);

		if (includedFields["DeletionQV"] and seq.deletionQV.Empty() == false) {
			deletionQVArray.Write(seq.deletionQV.data, seq.length);
		}
		if (includedFields["PreBaseDeletionQV"] and seq.preBaseDeletionQV.Empty() == false) {
			DNALength readPos;
			for (readPos = 0; readPos < seq.length; readPos++) {
				preBaseDeletionQVArray.WriteRow(&seq.preBaseDeletionQV[readPos*4], 4);
			}
		}
		if (includedFields["DeletionTag"] and seq.deletionTag != NULL) {
			deletionTagArray.Write(seq.deletionTag, seq.length);
		}
		if (includedFields["InsertionQV"] and seq.insertionQV.Empty() == false) {
			insertionQVArray.Write(seq.insertionQV.data, seq.length);
		}
		if (includedFields["SubstitutionQV"] and seq.substitutionQV.Empty() == false) {
			substitutionQVArray.Write(seq.substitutionQV.data, seq.length);
		}
		if (includedFields["SubstitutionTag"] and seq.substitutionTag != NULL) {
			substitutionTagArray.Write(seq.substitutionTag, seq.length);
		}
		if (includedFields["MergeQV"] and seq.mergeQV.Empty() == false) {
			mergeQVArray.Write(seq.mergeQV.data, seq.length);
		}

	}

	int WriteBases(FASTASequence &seq ) {
		int lenArray[1] = {seq.length};
		nElemArray.Write(lenArray, 1);
		baseArray.Write((const unsigned char*) seq.seq, seq.length);
		return 1;
	}

	int Write(SMRTSequence &seq) {
		WriteBases(seq);
		WriteQualities(seq);

		if (includedFields["PreBaseFrames"] and seq.preBaseFrames != NULL) {
			preBaseFramesArray.Write(seq.preBaseFrames, seq.length);
		}
		if (includedFields["WidthInFrames"] and seq.widthInFrames != NULL) {
			widthInFramesArray.Write(seq.widthInFrames, seq.length);
		}
    if (includedFields["PulseIndex"] and seq.pulseIndex != NULL) {
      pulseIndexArray.Write(seq.pulseIndex, seq.length);
    }
		WriteIdentifiers(seq.zmwData.holeNumber, seq.zmwData.holeStatus,  seq.xy[0], seq.xy[1]);

		return 1;
	}

	int Write(FASTQSequence &seq) {

		int x, y;
		UInt holeNumber;
    unsigned char holeStatus;

		WriteBases(seq);
		WriteQualities(seq);

		if (platformId == Astro) {
			//
			// now extract the x an y coordinates.
			GetSMRTReadCoordinates(seq, x, y);
      holeStatus = 0;
			seq.GetHoleNumber((int&) holeNumber);
		}
		if( platformId == Springfield){ 
			GetSpringfieldHoleNumberFromTitle(seq, holeNumber);
		}

		WriteIdentifiers(holeNumber,holeStatus,x,y);

		return 1;
	}
};





#endif
