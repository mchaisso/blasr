#ifndef DATA_HDF_HDF_PULSE_DATA_FILE_H_
#define DATA_HDF_HDF_PULSE_DATA_FILE_H_

#include "HDFGroup.h"
#include "HDFScanDataReader.h"

class HDFPulseDataFile {
 public:
	H5File hdfBasFile;
	HDFGroup pulseDataGroup;
	HDFGroup rootGroup;
	HDFGroup *rootGroupPtr;
	string pulseDataGroupName;
	HDFScanDataReader scanDataReader;
	bool useScanData;
	bool closeFileOnExit;
	int  maxAllocNElements;
	HDFZMWReader zmwReader;
  vector<unsigned int> eventOffset;
	int nReads;
  bool preparedForRandomAccess;
  map<unsigned int,int> holeNumberToIndex;
	vector<unsigned int> holeNumbers;
	vector<unsigned int> readLengths;

	int GetAllReadLengths(vector<DNALength> &readLengths) {
		nReads = zmwReader.numEventArray.arrayLength;
		readLengths.resize(nReads);
		zmwReader.numEventArray.Read(0,nReads, (int*) &readLengths[0]);
		return readLengths.size();
	}

	void CheckMemoryAllocation(long allocSize, long allocLimit, const char *fieldName = NULL) {
		if (allocSize > allocLimit) {
			if (fieldName == NULL) {
				cout << "Allocating too large of memory" << endl;
			}
			else {
                cout << "Allocate size " << allocSize << " > allocate limit " << allocLimit << endl;
				cout << "ERROR! Reading the dataset " << fieldName << " will use too much memory." << endl;
				cout << "The pls/bas file is too large, exiting." << endl;
			}
			exit(1);
		}
	}
	
	HDFPulseDataFile() {
		pulseDataGroupName = "PulseData";
		nReads             = 0;
		useScanData        = false;
		closeFileOnExit    = false;
		maxAllocNElements  = INT_MAX;
        preparedForRandomAccess = false;
        rootGroupPtr       = NULL;
	}

  void PrepareForRandomAccess() {
    GetAllReadLengths(readLengths);
		GetAllHoleNumbers(holeNumbers);
    int i;
    int curOffset = 0;
		eventOffset.resize(readLengths.size());
    for (i = 0; i < readLengths.size(); i++) {
      eventOffset[i] = curOffset;
      curOffset = curOffset + readLengths[i];
    }
		for (i = 0; i < holeNumbers.size(); i++) {
			holeNumberToIndex[holeNumbers[i]] = i;
		}
    nReads = eventOffset.size();
    preparedForRandomAccess = true;
  }


	int OpenHDFFile(string fileName, const H5::FileAccPropList & fileAccPropList=H5::FileAccPropList::DEFAULT) {
		try {
			H5::FileAccPropList propList = fileAccPropList;
      Exception::dontPrint();
			hdfBasFile.openFile(fileName.c_str(), H5F_ACC_RDONLY, propList);	
		}
		catch (Exception &e) {
			cout << "ERROR, could not open hdf file" << fileName << ", exiting." << endl;
			exit(1);
		}
		closeFileOnExit = true;
		return 1;
	}

	//
	// All pulse data files contain the "PulseData" group name.
	// 
	//
	int InitializePulseDataFile(string fileName, const H5::FileAccPropList & fileAccPropList=H5::FileAccPropList::DEFAULT) {
		if (OpenHDFFile(fileName, fileAccPropList) == 0) return 0;
		return 1;
	}

	int Initialize(string fileName, const H5::FileAccPropList & fileAccPropList=H5::FileAccPropList::DEFAULT) {

		if (InitializePulseDataFile(fileName, fileAccPropList) == 0) {
			return 0;
		}
		//
		// The pulse group is contained directly below the root group.
		//
		if (rootGroup.Initialize(hdfBasFile, "/") == 0) {
			return 0;
		}
		rootGroupPtr = &rootGroup;
		return Initialize();
	}

	//
	// Initialize inside another open group.
	//
	int Initialize(HDFGroup *rootGroupP) {
		rootGroupPtr = rootGroupP;
		return Initialize();
	}

	//
	// Initialize all fields 
    //
	int Initialize() {
    preparedForRandomAccess = false;		
		if (InitializePulseGroup() == 0) return 0;
		if (rootGroupPtr->ContainsObject("ScanData") == false or
				scanDataReader.Initialize(rootGroupPtr) == 0) {
			useScanData = false;
		}
		else {
			useScanData = true;
		}
		return 1;
	}

	int InitializePulseGroup() {
		if (pulseDataGroup.Initialize(rootGroupPtr->group, pulseDataGroupName) == 0) return 0;
		return 1;
	}

	int GetAllHoleNumbers(vector<unsigned int> &holeNumbers) {
		CheckMemoryAllocation(zmwReader.holeNumberArray.arrayLength, maxAllocNElements, "HoleNumbers (base)");
		holeNumbers.resize(nReads);
		zmwReader.holeNumberArray.Read(0,nReads, (unsigned int*)&holeNumbers[0]);
		return holeNumbers.size();
	}	

	void Close() {
		if (useScanData) {
			scanDataReader.Close();
		}
		
		pulseDataGroup.Close();
		if (rootGroupPtr == &rootGroup) {
			rootGroup.Close();
		}
		/*
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_FILE) << " open files upon closing." <<endl;
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_DATASET) << " open datasets upon closing." <<endl;
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_GROUP) << " open groups upon closing." <<endl;
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_DATATYPE) << " open datatypes upon closing." <<endl;
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_ATTR) << " open attributes upon closing." <<endl;
		*/
		if (closeFileOnExit) {
			hdfBasFile.close();
		}


	}

};

#endif
