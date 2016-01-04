#ifndef DATA_HDF_HDF_BAS_READER_H_
#define DATA_HDF_HDF_BAS_READER_H_

#include <stdlib.h>
#include <sstream>
#include <vector>

#include "DatasetCollection.h"
#include "HDFArray.h"
#include "HDF2DArray.h"
#include "HDFAtom.h"
#include "HDFGroup.h"
#include "HDFZMWReader.h"
#include "HDFScanDataReader.h"
#include "HDFPulseDataFile.h"
#include "datastructures/reads/BaseFile.h"
#include "utils/VectorUtils.h"
#include "FASTQSequence.h"
#include "SMRTSequence.h"
#include "Enumerations.h"
#include "utils/ChangeListID.h"
using namespace H5;
using namespace std;

/*
 * Below is sample code for using the bas reader to read in sequences
 * from a .bas.h5 or .pls.h5 file.  
 * One may select which quality value fields to read.  By default,
 * none are read in.  To read in the rich quality values (insertionQV,
 * deletionQV, substitutionQV, mergeQV, substitutionTag, deletionTag),
 * call InitializeDefaultIncludedFields() BEFORE initializing the
 * reader on the file name.
 *
#include "data/hdf/HDFBasReader.h"
#include "SMRTSequence.h"
#include <string>

int main(int argc, char* argv[]) {
	
	string basFile = argv[1];

	HDFBasReader reader;
  reader.InitializeDefaultIncludedFields();
	reader.Initialize(basFile);

	SMRTSequence read;

	while(reader.GetNext(read)) {
		read.PrintQualSeq(cout);
	}


	return 0;
}

 * If you wanted to read only fasta sequences and not fastq, use:
 
 FASTASequence read;
 while (reader.GetNext(read) {
   read.PritnSeq(cout);
 }

 */

template<typename T_Sequence>
class T_HDFBasReader : public DatasetCollection, public HDFPulseDataFile {
 public:
	DNALength curBasePos;
	int curRead;
	unsigned int nBases;
  

	//bool readPulseInformation;
	bool hasRegionTable;

	HDFArray<int> zmwXCoordArray;
	HDFArray<int> zmwYCoordArray;
	HDFArray<unsigned char> baseArray;
	HDFArray<unsigned char> deletionQVArray;
	HDFArray<unsigned char> deletionTagArray;
	HDFArray<unsigned char> insertionQVArray;
	HDFArray<unsigned char> substitutionTagArray;
	HDFArray<unsigned char> substitutionQVArray;
	HDFArray<unsigned char> mergeQVArray;
  HDFArray<unsigned char> qualArray;
	HDFArray<unsigned int> simulatedCoordinateArray;
	HDFArray<unsigned int> simulatedSequenceIndexArray;
	HDFArray<uint16_t> basWidthInFramesArray;
	HDFArray<uint16_t> preBaseFramesArray;
	HDFArray<int> pulseIndexArray;
	HDFArray<int> holeIndexArray;
	HDFGroup baseCallsGroup;
	HDFGroup zmwGroup;
//	HDFArray<HalfWord> pulseWidthArray; This is deprecated
  HDFAtom<string> changeListIDAtom;

	//bool useWidthInFrames,  usePulseIndex, 
    bool useZmwReader;
//    bool usePulseWidth; This is deprecated
	
	string baseCallsGroupName;
	bool qualityFieldsAreCritical;
	
	//bool useHoleNumbers, useHoleStatus;
	//bool usePreBaseFrames;
	bool useBasHoleXY;
	//bool useBasecall;
	//bool useQuality;
  bool readBasesFromCCS;
  QVScale qvScale;

 public:
	PlatformId GetPlatform() {
		return scanDataReader.platformId;
	}

	string GetMovieName() {
		if (scanDataReader.useMovieName) {
			return scanDataReader.GetMovieName();
		}
		else {
			return "";
		}
	}

  int GetReadAt(int holeNumber, SMRTSequence &read) {
    //
    // The first time this is called there may be no table of read
    // offset positions.  Check for that and build if it does not
    // exist. 
    //
    if (preparedForRandomAccess == false) {
      PrepareForRandomAccess();
    }
		int index  = holeNumberToIndex[holeNumber];
    curRead    = index;
    curBasePos = eventOffset[index];
    zmwReader.curZMW = index;
    return GetNext(read);
  }

	
	string GetRunCode() {
		return scanDataReader.GetRunCode();
	}

	T_HDFBasReader() {
		curRead      = 0;
		curBasePos   = 0;
        nBases       = 0;
        preparedForRandomAccess = false;
        readBasesFromCCS = false;
		baseCallsGroupName = "BaseCalls";
		qualityFieldsAreCritical = true;
		useZmwReader = false;
        useBasHoleXY = true;
        hasRegionTable = false;
        qvScale = POverOneMinusP; //default 0 = POverOneMinusP
		fieldNames.push_back("Basecall");
		fieldNames.push_back("DeletionQV");
		fieldNames.push_back("DeletionTag");
		fieldNames.push_back("InsertionQV");
		fieldNames.push_back("SubstitutionTag");
		fieldNames.push_back("SubstitutionQV");
		fieldNames.push_back("QualityValue");
		fieldNames.push_back("WidthInFrames");
		fieldNames.push_back("PulseIndex");
		fieldNames.push_back("PreBaseFrames");
		fieldNames.push_back("MergeQV");
		fieldNames.push_back("SimulatedCoordinate");
		fieldNames.push_back("SimulatedSequenceIndex");
    // Start out with no fields being read.
		InitializeAllFields(false);
    // Then by default always read bases.
        IncludeField("Basecall");
	}


  void InitializeDefaultCCSIncludeFields() {
		InitializeAllFields(false);
		IncludeField("Basecall");
		IncludeField("DeletionQV");
		IncludeField("DeletionTag");
		IncludeField("InsertionQV");
		IncludeField("SubstitutionQV");
		IncludeField("SubstitutionTag");
		IncludeField("QualityValue");
  }

  void InitializeDefaultRawBasIncludeFields() {
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
    IncludeField("MergeQV");
  }
	void InitializeDefaultIncludedFields() {
    if (readBasesFromCCS == false) {
      InitializeDefaultRawBasIncludeFields();
    }
    else { 
      InitializeDefaultCCSIncludeFields();
    }
	}
	
	void InitializeDefaultRequiredFields() {
		requiredFields["Basecall"] = true;
	}

	bool HasRegionTable() {
		return hasRegionTable;
	}
	
	
	bool HasSimulatedCoordinatesStored() {
		if (baseCallsGroup.ContainsObject("SimulatedCoordinate") and
				baseCallsGroup.ContainsObject("SimulatedSequenceIndex")) {
			return true;
		}
		else {
			return false;
		}
	}
  
  void SetReadBasesFromCCS() {
    InitializeDefaultCCSIncludeFields();
    readBasesFromCCS = true;
  }
  
  void GetChangeListID(string &changeListID) {
    if (changeListIDAtom.initialized) {
      changeListIDAtom.Read(changeListID);
    }
    else {
      changeListID = "0";
    }
  }

	void GetChangelistId(string &changelistId) {
      changeListIDAtom.Initialize(baseCallsGroup.group, "ChangeListID");
      GetChangeListID(changelistId);
	}

	int InitializeForReadingBases() {

		//
		// Initialize root group + scan data information.
		//
		if (HDFPulseDataFile::Initialize(rootGroupPtr) == 0) return 0;
    
		//
		// Open the base group, this contains all the required information.  
		//

    if (readBasesFromCCS) {
      baseCallsGroupName = "ConsensusBaseCalls";
    }
    if (pulseDataGroup.ContainsObject(baseCallsGroupName) == 0 or
        baseCallsGroup.Initialize(pulseDataGroup.group, baseCallsGroupName) == 0) {
      return 0;
    }
    if (baseCallsGroup.ContainsAttribute("ChangeListID")) {
      changeListIDAtom.Initialize(baseCallsGroup.group, "ChangeListID");
      string changeListIdString;
      ChangeListID changeList;
      GetChangeListID(changeListIdString);
      qvScale = DetermineQVScaleFromChangeListID(changeList);
    }
		if (pulseDataGroup.ContainsObject("Regions")) {
			hasRegionTable = true;
		}
		else {
			hasRegionTable = false;
		}

		//
		// Initialize read and quality arrays for reading.
		//
		if (this->InitializeSequenceFields(baseCallsGroup) == 0) {
      return 0;
    }

		//
		// Initialize simulated coordinate fields if they exist.  They are
		// automatically opened and read from when they are in the bas.h5
		// file since they are used for debugging.
		//
		
		if (baseCallsGroup.ContainsObject("SimulatedCoordinate")) {
			includedFields["SimulatedCoordinate"] = true;
			InitializeDataset(baseCallsGroup, simulatedCoordinateArray, "SimulatedCoordinate");
		}
		else {
			includedFields["SimulatedCoordinate"] = false;
		}

		if (baseCallsGroup.ContainsObject("SimulatedSequenceIndex")) {
			includedFields["SimulatedSequenceIndex"] = true;
			InitializeDataset(baseCallsGroup, simulatedSequenceIndexArray, "SimulatedSequenceIndex");
		}
		else {
			includedFields["SimulatedSequenceIndex"] = false;
		}
		nBases = baseArray.arrayLength;

		return 1;
	}


	int InitializeCommon() {

		//
		// Initialize the smallest set of fields required to import bases.
		//
		if (InitializeForReadingBases() == 0) {
			return 0;
		}
		

		return 1;
	}

	template<typename T_Dataset>
	void InitializeRequiredField(HDFGroup &group, string arrayName, T_Dataset &field) {
		if (group.ContainsObject(arrayName)) {
			if (field.Initialize(group, arrayName) != 0) {
				return;
			}
		}
		cout << "ERROR. Could not initialize dataset " << arrayName << endl;
		exit(1);
	}
		
	template<typename T>
	int InitializeField(HDFGroup &rootGroup, string arrayName,
											T &field, bool &initialized) {
		initialized = false;
		if (rootGroup.ContainsObject(arrayName)) {
			if (field.Initialize(rootGroup, arrayName) != 0) {
				initialized = true;
				return true;
			}
		}
		return false;
	}

	template<typename T>
	int InitializeAttribute(HDFGroup &rootGroup, string arrayName,
											T &field, bool &initialized,
											bool fieldIsCritical = true) {
		int success = 1;
		initialized = false;
		if (rootGroup.ContainsAttribute(arrayName)) {
			if (field.Initialize(rootGroup, arrayName) == 0) {
				success = 0;
			}
			else {
				initialized = true;
			}
		}
		else {
			// the field does not exist
			success = 0;
		}
		if (fieldIsCritical) {
			return success;
		}
		else {
			return 1;
		}
	}
														 
	int InitializeSequenceFields(HDFGroup &baseCallsGroup) {
		//
		// The only field that is absoultely required is Basecall
		if (baseArray.InitializeForReading(baseCallsGroup, "Basecall")        == false) return 0;
		if (includedFields["QualityValue"] and qualArray.InitializeForReading(baseCallsGroup, "QualityValue")    == false) return 0;
		if (includedFields["InsertionQV"] and insertionQVArray.InitializeForReading(baseCallsGroup, "InsertionQV")     == false) return 0;
		if (includedFields["DeletionQV"] and deletionQVArray.InitializeForReading(baseCallsGroup, "DeletionQV")      == false) return 0;
		if (includedFields["DeletionTag"] and deletionTagArray.InitializeForReading(baseCallsGroup, "DeletionTag")     == false) return 0;
		if (includedFields["SubstitutionQV"] and substitutionQVArray.InitializeForReading(baseCallsGroup, "SubstitutionQV")  == false) return 0;
		if (includedFields["SubstitutionTag"] and substitutionTagArray.InitializeForReading(baseCallsGroup, "SubstitutionTag") == false) return 0;
		// if (includedFields["PreBaseFrames"] and preBaseFramesArray.InitializeForReading(baseCallsGroup, "PreBaseFrames")   == false) return 0;
        
        if (baseCallsGroup.ContainsObject("PreBaseFrames")) {
            if (preBaseFramesArray.InitializeForReading(baseCallsGroup, "PreBaseFrames") == false) return 0;
        } else {
            includedFields["PreBaseFrames"] = false;
        }
    

    //
    // These fields are not always present in bas.h5 files.
    //
    if (baseCallsGroup.ContainsObject("PulseIndex")) {
      if (pulseIndexArray.InitializeForReading(baseCallsGroup,        "PulseIndex")      == false) return 0;
    }
    else {
      includedFields["PulseIndex"] = false;
    }

    if (baseCallsGroup.ContainsObject("WidthInFrames")) {    
      if (basWidthInFramesArray.InitializeForReading(baseCallsGroup,  "WidthInFrames")   == false) return 0;
    }
    else {
      includedFields["WidthInFrames"] = false;
    }

    if (baseCallsGroup.ContainsObject("MergeQV")) {
      if (includedFields["MergeQV"] and mergeQVArray.InitializeForReading(baseCallsGroup, "MergeQV") == false) return false;
    }
    else {
      includedFields["MergeQV"] = false;
    }


    return 1;
	}

  int InitializeAstro() {
		useBasHoleXY = true;
		return 1;
	}

	int InitializeSpringfield() {
		//
		// For now, no special initialization is required.
		//
		return 1;
	}

	int Initialize(HDFGroup *rootGroupP) {
		rootGroupPtr= rootGroupP;
		return Initialize();
	}


	int Initialize() {
		
		// Return 0 if any of the array inializations do not work.
		if (InitializeCommon() == 0) {
			return 0;
		}

		// Must have zmw information.
		if (zmwReader.Initialize(&baseCallsGroup) == 0) {
			return 0;
		}
		else {
			useZmwReader = true;
		}
		//
		// Get information about the chip - how many zmw's, the hole
		// number indices to keep track of reads, etc..
		//
		nReads = zmwReader.numEventArray.arrayLength;
		
		

		if (scanDataReader.platformId == AstroPlatform) {
			if (InitializeAstro() == 0) {
				return 0;
			}
		}
		else if (scanDataReader.platformId == SpringfieldPlatform) {
			if (InitializeSpringfield() == 0) {
				return 0;
			}
		}

		/*
		 * Initialize state variables.
		 */
		
		curBasePos = 0;
		curRead    = 0;

		/*
		 * All ok, return that.
		 */
		return 1;
	}

    int InitializeHDFFile(string hdfBasFileName, 
            const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT) { 
		/*
		 * Initialize access to the HDF file.  For reading bas files, this
		 * involves:
		 *   - Opening the file, and initializing both base and zmw grops.
		 */
		if (OpenHDFFile(hdfBasFileName, fileAccPropList) == 0) {
            return 0;
		}

		if (rootGroup.Initialize(hdfBasFile, "/") == 0) {
			return 0;
		}
		rootGroupPtr = &rootGroup;
        return 1;
    }

	int Initialize(string hdfBasFileName,
            const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT) {
        int init = InitializeHDFFile(hdfBasFileName, fileAccPropList);
        if (init == 0) return 0;
		return Initialize();
    }
   
	int GetNumReads() {
		return nReads;
	}

	void BuildReadTitle(string movieTitle, unsigned int holeNumber, string &readTitle, unsigned int simIndex=0, unsigned int simCoordinate=0) {
		stringstream readTitleStrm;
		readTitleStrm << movieTitle << "/" << holeNumber;
		readTitle = readTitleStrm.str();
	}

	int GetNext(FASTASequence &seq) {
		if (curRead == nReads) {
			return 0;
		}
	
		int seqLength;
        try {
		seqLength = GetNextWithoutPosAdvance(seq);
        } catch(DataSetIException e) {
            cout << "ERROR, could not read base calls for FASTA Sequence "
                 << seq.GetName() << endl;
            exit(1);
        }
		curBasePos += seqLength;
		seq.StorePlatformType(scanDataReader.platformId);
		return 1;
	}

	
	int GetNext(FASTQSequence &seq) {
        try {
		if (curRead == nReads) {
			return 0;
		}
		int seqLength = GetNextWithoutPosAdvance(seq);
		seq.length = seqLength;

		if (seqLength > 0 ) {
			if (includedFields["QualityValue"]) {
				seq.AllocateQualitySpace(seqLength);
				qualArray.Read((int)curBasePos, (int) curBasePos + seqLength, (unsigned char*) seq.qual.data);
			}
		}

		if (includedFields["DeletionQV"]) {
			GetNextDeletionQV(seq);
		}
		if (includedFields["DeletionTag"]) {
			GetNextDeletionTag(seq);
		}
		if (includedFields["InsertionQV"]) {
			GetNextInsertionQV(seq);
		}
		if (includedFields["SubstitutionQV"]) {
			GetNextSubstitutionQV(seq);
		}
		if (includedFields["SubstitutionTag"]) {
			GetNextSubstitutionTag(seq);
        }
        if (includedFields["MergeQV"]) {
            GetNextMergeQV(seq);
        }
        seq.SetQVScale(qvScale);
		curBasePos += seqLength;
        } catch(DataSetIException e) {
            cout << "ERROR, could not read quality metrics for FASTQ Sequence " 
                 << seq.GetName() << endl;
            exit(1);
        }
        return 1;
	}

//
// Reading of SMRT Sequences reads both the sequence fields, and the
// fields with ZMW information for identification of this read.
//

 int GetNext(SMRTSequence &seq) {
	 //
	 // Read in quality values.
	 //
	 int retVal;
	 
	 DNALength  curBasPosCopy = curBasePos;
	 //
	 // Getting next advances the curBasPos to the end of 
	 // the current sequence. 
	 //

	 retVal = this->GetNext((FASTQSequence&)seq);
	 //
	 // Bail now if the file is already done
     //
	 if (retVal  == 0) {
		 return 0;
	 }

     try {
	 DNALength nextBasePos = curBasePos;
	 curBasePos = curBasPosCopy;

	 if (includedFields["WidthInFrames"] ) {
		 assert(nextBasePos <= basWidthInFramesArray.arrayLength);
		 GetNextWidthInFrames(seq);
	 }
	 if (includedFields["PreBaseFrames"]) { 
		 GetNextPreBaseFrames(seq);
	 }
	 if (includedFields["PulseIndex"]) { 
		 GetNextPulseIndex(seq);
	 }
     	 curBasePos = nextBasePos;
	 
	 //
	 // By default, the subread of a read without subread information is
	 // the whole read.
	 //
	 seq.subreadStart = 0;
	 seq.subreadEnd   = seq.length;
     zmwReader.GetNext(seq.zmwData);
     seq.xy[0] = seq.zmwData.x;
     seq.xy[1] = seq.zmwData.y;
     } catch(DataSetIException e) {
         cout << "ERROR, could not read pulse metrics for SMRTSequence " 
              << seq.GetName() << endl;
         exit(1);
     }
	 return retVal;
 }
 /*
    int16_t xy[2];
    if (zmwReader.readHoleXY) {
      zmwReader.xyArray.Read(curRead, curRead+1, 0, 2, xy);
    }
    else {
      xy[0] = xy[1] = 0;
    }
    seq.StoreXY(xy);
*/
	void GetAllPulseIndex(vector<int> &pulseIndex) {
		CheckMemoryAllocation(pulseIndexArray.arrayLength, maxAllocNElements, "PulseIndex");
		pulseIndex.resize(pulseIndexArray.arrayLength);
		pulseIndexArray.Read(0, pulseIndexArray.arrayLength, &pulseIndex[0]);
	}

	int GetAllPreBaseFrames(vector<uint16_t> &preBaseFrames) {
		CheckMemoryAllocation(preBaseFramesArray.arrayLength, maxAllocNElements, "PreBaseFrames");
		preBaseFrames.resize(nBases);
		preBaseFramesArray.Read(0, nBases, &preBaseFrames[0]);
	}

	int GetAllWidthInFrames(vector<uint16_t> &widthInFrames) { 
		CheckMemoryAllocation(basWidthInFramesArray.arrayLength, maxAllocNElements, "WidthInFrames");
		widthInFrames.resize(nBases);
		basWidthInFramesArray.Read(0, nBases, &widthInFrames[0]);
	}

  int GetAllHoleStatus(vector<unsigned char> &holeStatus) {
		CheckMemoryAllocation(zmwReader.holeStatusArray.arrayLength, maxAllocNElements, "HoleStatus (base)");
		holeStatus.resize(nReads);
		zmwReader.holeStatusArray.Read(0,nReads, (unsigned char*)&holeStatus[0]);
		return holeStatus.size();
	}

  int GetAllReadLengths(vector<int> &readLengths) {
    readLengths.resize(nReads);
		zmwReader.numEventArray.ReadDataset(readLengths);
    return readLengths.size();
  }

	int Advance(int nSeq) {
		int i;
		// cannot advance past the end of this file
		if (curRead + nSeq >= nReads) { return 0; }
		for (i = curRead; i < curRead + nSeq && i < nReads; i++ ) {
 			int seqLength;
			zmwReader.numEventArray.Read(i, i+1, &seqLength);
			curBasePos += seqLength;
		}
		curRead += nSeq;
		zmwReader.Advance(nSeq);
		return curRead;
	}

	int GetNextWithoutPosAdvance(FASTASequence &seq) {
		int seqLength;

		zmwReader.numEventArray.Read(curRead, curRead+1, &seqLength);
		seq.length = 0;
		seq.seq = NULL;

		if (includedFields["Basecall"]) {
			if (seqLength > 0) {
				ResizeSequence(seq, seqLength);
				baseArray.Read(curBasePos, curBasePos + seqLength, (unsigned char*) seq.seq);
			}
		}

		string readTitle;
		unsigned int holeNumber;
    unsigned char holeStatus;
		zmwReader.holeNumberArray.Read(curRead, curRead+1, &holeNumber);
		seq.StoreHoleNumber(holeNumber);
		seq.StoreHoleStatus(holeStatus);

		DNALength simIndex=0, simCoordinate=0;

		if (includedFields["SimulatedSequenceIndex"] == true) {
			simulatedSequenceIndexArray.Read(curRead,curRead+1,&simIndex);
		}
		if (includedFields["SimulatedCoordinate"] == true) {
			simulatedCoordinateArray.Read(curRead, curRead+1, &simCoordinate);
		}
		


    BuildReadTitle(scanDataReader.GetMovieName(), holeNumber, readTitle, simIndex, simCoordinate);

		seq.CopyTitle(readTitle);
		curRead++;
		return seqLength;
	}

	int GetNextDeletionQV(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateDeletionQVSpace(seq.length);
		deletionQVArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.deletionQV.data);
    return seq.length;
	}

	int GetNextMergeQV(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateMergeQVSpace(seq.length);
		mergeQVArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.mergeQV.data);
    return seq.length;
	}

	int GetNextDeletionTag(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateDeletionTagSpace(seq.length);
		deletionTagArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.deletionTag);
    return seq.length;
	}

	int GetNextInsertionQV(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateInsertionQVSpace(seq.length);
		insertionQVArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.insertionQV.data);
    return seq.length;
	}

	int GetNextWidthInFrames(SMRTSequence &seq) {
		if (seq.length == 0) return 0;
		seq.widthInFrames = new HalfWord[seq.length];
		basWidthInFramesArray.Read((int)curBasePos, (int) curBasePos + seq.length, (HalfWord*) seq.widthInFrames);
    return seq.length;
	}

	int GetNextPreBaseFrames(SMRTSequence &seq) {
		if (seq.length == 0) return 0;
		seq.preBaseFrames = new HalfWord[seq.length];
		preBaseFramesArray.Read((int)curBasePos, (int) curBasePos + seq.length, (HalfWord*) seq.preBaseFrames);
    return seq.length;
	}
	int GetNextPulseIndex(SMRTSequence &seq) {
		if (seq.length == 0) return 0;
		seq.pulseIndex = new int[seq.length];
		pulseIndexArray.Read((int)curBasePos, (int) curBasePos + seq.length, (int*) seq.pulseIndex);
    return seq.length;
	}

	int GetNextSubstitutionQV(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateSubstitutionQVSpace(seq.length);
		substitutionQVArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.substitutionQV.data);
    return seq.length;
	}

	int GetNextSubstitutionTag(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateSubstitutionTagSpace(seq.length);
		substitutionTagArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.substitutionTag);		
    return seq.length;
	}

	void Close() {

		baseCallsGroup.Close();
		zmwXCoordArray.Close();
		zmwYCoordArray.Close();
		baseArray.Close();
		qualArray.Close();
		if (useZmwReader) {
			zmwReader.Close();
		}

		if (includedFields["DeletionQV"]) {
			deletionQVArray.Close();
		}
		if (includedFields["DeletionTag"]) {
			deletionTagArray.Close();
		}
		if (includedFields["MergeQV"]) {
			mergeQVArray.Close();
		}
		if (includedFields["InsertionQV"]) {
			insertionQVArray.Close();
		}
		if (includedFields["SubstitutionTag"]) {
			substitutionTagArray.Close();
		}
		if (includedFields["SubstitutionQV"]) {
			substitutionQVArray.Close();
		}
		if (includedFields["WidthInFrames"]) {
			basWidthInFramesArray.Close();
		}
		if (includedFields["PreBaseFrames"]) {
			preBaseFramesArray.Close();
		}
		if (includedFields["PulseIndex"]) {
			pulseIndexArray.Close();
		}

		HDFPulseDataFile::Close();
	}

	void ReadAllHoleXY(BaseFile &baseFile) {
		baseFile.holeXY.resize(nReads);
		int i;
		for (i = 0; i < nReads; i++) {
			zmwReader.xyArray.Read(i,i+1, baseFile.holeXY[i].xy);
		}
	}

    //
    // Return size of an entire field.
    //
    UInt GetFieldSize(const string & field) {
        if (not includedFields[field]) {
            cout << "ERROR, field [" << field << "] is not included in the base file." << endl;
            exit(1);
        }
        if (field == "Basecall") {
            return baseArray.arrayLength / 1024 * sizeof(unsigned char);
        } else if (field == "QualityValue") {
            return qualArray.arrayLength / 1024 * sizeof(unsigned char);
        } else if (field == "DeletionQV") {
            return deletionQVArray.arrayLength / 1024 * sizeof(unsigned char);
        } else if (field == "DeletionTag")  {
            return deletionTagArray.arrayLength / 1024 * sizeof(unsigned char);
        } else if (field == "MergeQV") {
            return mergeQVArray.arrayLength / 1024 * sizeof(unsigned char);
        } else if (field == "InsertionQV") { 
            return insertionQVArray.arrayLength / 1024 * sizeof(unsigned char);
        } else if (field == "SubstitutionQV") {
            return substitutionQVArray.arrayLength / 1024 * sizeof(unsigned char);
        } else if (field == "SubstitutionTag") {
            return substitutionTagArray.arrayLength / 1024 * sizeof(unsigned char);
        } else if (field == "WidthInFrames") {
            return basWidthInFramesArray.arrayLength / 1024 * sizeof(uint16_t);
        } else if (field == "PreBaseFrames") {
			return preBaseFramesArray.arrayLength / 1024 * sizeof(uint16_t);
        } else if (field == "PulseIndex") {
            return pulseIndexArray.arrayLength / 1024 * sizeof(int);
        } else {
            cout << "ERROR, field [" << field << "] is not supported. " << endl ;
            exit(1);
        }
    }

    //
    // Read an entire field.
    //
    void ReadField(BaseFile & baseFile, const string & field) {
        if (not includedFields[field]) {
            cout << "ERROR, field [" << field << "] is not included in the base file." << endl;
            exit(1);
        }
        if (field == "Basecall") {
            assert(nBases == baseArray.arrayLength);
            baseArray.ReadDataset(baseFile.baseCalls); 
        } else if (field == "QualityValue") {
            qualArray.ReadDataset(baseFile.qualityValues);
        } else if (field == "DeletionQV") {
            deletionQVArray.ReadDataset(baseFile.deletionQV);
        } else if (field == "DeletionTag")  {
            deletionTagArray.ReadDataset(baseFile.deletionTag);
        } else if (field == "MergeQV") {
            mergeQVArray.ReadDataset(baseFile.mergeQV);
        } else if (field == "InsertionQV") { 
			insertionQVArray.ReadDataset(baseFile.insertionQV);
        } else if (field == "SubstitutionQV") {
			substitutionQVArray.ReadDataset(baseFile.substitutionQV);			
        } else if (field == "SubstitutionTag") {
            substitutionTagArray.ReadDataset(baseFile.substitutionTag);
        } else if (field == "WidthInFrames") {
			basWidthInFramesArray.ReadDataset(baseFile.basWidthInFrames);
        } else if (field == "PreBaseFrames") {
			preBaseFramesArray.ReadDataset(baseFile.preBaseFrames);
        } else if (field == "PulseIndex") {
			pulseIndexArray.ReadDataset(baseFile.pulseIndex);
        } else {
            cout << "ERROR, field [" << field << "] is not supported. " << endl ;
            exit(1);
        }
    }

    //
    // Clear memory allocated for a field
    //
    void ClearField(BaseFile & baseFile, const string & field) {
        if (not includedFields[field]) {
            cout << "ERROR, field [" << field << "] is not included in the base file." << endl;
            exit(1);
        }
        if (field == "Basecall") {
    	    ClearMemory(baseFile.baseCalls);
        } else if (field == "QualityValue") {
            ClearMemory(baseFile.qualityValues);
        } else if (field == "DeletionQV") {
            ClearMemory(baseFile.deletionQV);
        } else if (field == "DeletionTag")  {
            ClearMemory(baseFile.deletionTag);
        } else if (field == "MergeQV") {
            ClearMemory(baseFile.mergeQV); 
        } else if (field == "InsertionQV") { 
            ClearMemory(baseFile.insertionQV);
        } else if (field == "SubstitutionQV") {
            ClearMemory(baseFile.substitutionQV);
        } else if (field == "SubstitutionTag") {
            ClearMemory(baseFile.substitutionTag);
        } else if (field == "WidthInFrames") {
            ClearMemory(baseFile.basWidthInFrames);
        } else if (field == "PreBaseFrames") {
            ClearMemory(baseFile.preBaseFrames); 
        } else if (field == "PulseIndex") {
            ClearMemory(baseFile.pulseIndex);
        } else {
            cout << "ERROR, field [" << field << "] is supported. " << endl ;
            exit(1);
        }
    }

    //
    // Initialization for reading a base file.  
    //
    void ReadBaseFileInit(BaseFile & baseFile) {
		if (scanDataReader.fileHasScanData) {
			scanDataReader.Read(baseFile.scanData);
		}

		baseFile.nReads = nReads;

		if (useBasHoleXY) {
			ReadAllHoleXY(baseFile);
		}
		GetAllHoleNumbers(baseFile.holeNumbers);
		GetAllHoleStatus(baseFile.holeStatus);
		zmwReader.numEventArray.ReadDataset(baseFile.readLengths);
		
        //
        // Cache the start positions of all reads.
        //
        assert(baseFile.nReads == baseFile.readLengths.size());
        baseFile.readStartPositions.resize(baseFile.readLengths.size()+1);

        if ( baseFile.readLengths.size() > 0 ) {
            int i;
            baseFile.readStartPositions[0] = 0;
            for (i = 1; i < baseFile.readLengths.size()+1; i++) {
                baseFile.readStartPositions[i] = 
                    (baseFile.readStartPositions[i-1] +
                     baseFile.readLengths[i-1]);
            }
        }
    }

    //
    // Read a base file.
    //
    void ReadBaseFile(BaseFile &baseFile) {
        ReadBaseFileInit(baseFile);

		if (includedFields["Basecall"]) {
			baseFile.baseCalls.resize(nBases);
			baseArray.Read(0,nBases, &baseFile.baseCalls[0]);
		}

		/*
		 * This can probably be fixed eventually with an object factory or
		 * collection of some sorts.
		 */
		if (includedFields["WidthInFrames"]) {
			basWidthInFramesArray.ReadDataset(baseFile.basWidthInFrames);
		}
		if (includedFields["PreBaseFrames"]) {
			preBaseFramesArray.ReadDataset(baseFile.preBaseFrames);
		}
		if (includedFields["PulseIndex"]) {
			pulseIndexArray.ReadDataset(baseFile.pulseIndex);
		}
		if (includedFields["QualityValue"]) {
			qualArray.ReadDataset(baseFile.qualityValues);
		}
		if (includedFields["InsertionQV"]) {
			insertionQVArray.ReadDataset(baseFile.insertionQV);
		}
		if (includedFields["SubstitutionTag"]) {
			substitutionTagArray.ReadDataset(baseFile.substitutionTag);
		}
		if (includedFields["SubstitutionQV"]) {
			substitutionQVArray.ReadDataset(baseFile.substitutionQV);			
		}
		if (includedFields["MergeQV"]) {
			mergeQVArray.ReadDataset(baseFile.mergeQV);			
		}
		if (includedFields["DeletionQV"]) {
			deletionQVArray.ReadDataset(baseFile.deletionQV);
		}
		if (includedFields["DeletionTag"]) {
			deletionTagArray.ReadDataset(baseFile.deletionTag);
		}

		baseFile.nBases = nBases;
		baseFile.scanData.platformId = scanDataReader.platformId;
	}
};


template<>
int T_HDFBasReader<SMRTSequence>::Advance(int nSteps) {
	int retVal;
	retVal = ((T_HDFBasReader<FASTQSequence>*)this)->Advance(nSteps);
	return retVal;
}

typedef T_HDFBasReader<FASTASequence> HDFBasReader;
typedef T_HDFBasReader<FASTQSequence> HDFQualReader;
typedef T_HDFBasReader<SMRTSequence>  HDFSmrtReader;


#endif
