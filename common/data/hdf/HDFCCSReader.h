#ifndef DATA_HDF_HDF_CCS_READER_H_
#define DATA_HDF_HDF_CCS_READER_H_

#include "HDFBasReader.h"

template<typename T_Sequence>
class HDFCCSReader : public T_HDFBasReader<T_Sequence> {
 public:
	HDFGroup ccsGroup, passesGroup;

	HDFArray<Byte> baseCallArray;
	HDFArray<UInt> passStartPulseArray, passNumPulsesArray, passStartBaseArray,
		passNumBasesArray, numPassesArray;
	HDFArray<Byte> passDirectionArray, adapterHitAfterArray, adapterHitBeforeArray;
	
	HDFZMWReader zmwReader;
	T_HDFBasReader<SMRTSequence> ccsBasReader;
	int curPassPos;

  HDFCCSReader() : T_HDFBasReader<T_Sequence>() {
		curPassPos = 0;
		this->fieldNames.push_back("AdapterHitAfter");
		this->fieldNames.push_back("AdapterHitBefore");
		this->fieldNames.push_back("NumPasses");
		this->fieldNames.push_back("PassDirection");
		this->fieldNames.push_back("PassNumPase");
		this->fieldNames.push_back("PassStartBase");
		this->fieldNames.push_back("PassStartPulse");
		this->fieldNames.push_back("PassNumPulses");
		InitializeAllCCSFields(true);
	}

	void InitializeAllCCSFields(bool value) {
		this->includedFields["AdapterHitAfter"]  = value;
		this->includedFields["AdapterHitBefore"] = value;
		this->includedFields["NumPasses"]        = value;
		this->includedFields["PassDirection"]    = value;
		this->includedFields["PassNumPase"]      = value;
		this->includedFields["PassStartBase"]    = value;
		this->includedFields["PassStartPulse"]   = value;
		this->includedFields["PassNumPulses"]    = value;
	}

	bool BasFileHasCCS(string ccsBasFileName) {
		try {
      Exception::dontPrint();
			this->hdfBasFile.openFile(ccsBasFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		}
		catch (H5::Exception &e) {
			cout << "ERROR, could not open hdf file " << ccsBasFileName << " Stopping." << endl;
			exit(1);
		}

		HDFGroup ccsBasecallsGroup;
		bool fileContainsCCS = false;
		HDFGroup pulseDataGroup;
		if (pulseDataGroup.Initialize(this->hdfBasFile, "PulseData") == 0) {
      cout << "ERROR, ccs base file " << ccsBasFileName 
           << " does not have a PulseData field." << endl;
      exit(1);
    }
		if (pulseDataGroup.ContainsObject("ConsensusBaseCalls")) {
			fileContainsCCS = true;
		}

		this->hdfBasFile.close();
		return fileContainsCCS;
	}

	int Advance(int nSteps) {
		cout << "ERROR! Advance is not yet implemented for ccs reader" << endl;
		assert(0);
		return 0;
	}

	int Initialize(string ccsBasFileName,
            const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT) {
		//
		// Open the file and initialize for reading reads.
		//

		// 
		// First, initialize for reading the unrolled bases from this
		// file.
		//

	  if (this->T_HDFBasReader<T_Sequence>::Initialize(ccsBasFileName, fileAccPropList) == 0) {
        cout << "ERROR, Could not initialize ccs file " << ccsBasFileName << endl;
				//     exit(1);
		}

		if (this->pulseDataGroup.ContainsObject("ConsensusBaseCalls")) {
			if (ccsGroup.Initialize(this->hdfBasFile, "PulseData/ConsensusBaseCalls") == 1) {
				/*
					cout << "ERROR, attempting to read cicular consensus data from '" << ccsBasFileName 
					<< "', which does not contain a ConsensusBaseCalls field." << endl;
					cout << "Check HDF file structure." << endl;
				*/
				//				exit(1);

				curPassPos = 0;
				int passesSuccess = 1;
				if (ccsGroup.ContainsObject("Passes") == 0) { 
					passesSuccess = 0;
				}
				else {
					if (passesGroup.Initialize(ccsGroup.group,"Passes") == 0) {
						passesSuccess = 0;
					}
				}
			
				if (passesSuccess == 0) {
					cout <<"ERROR, attempting to read circular consensus group Passes but it does not exist. " << endl;
					cout <<"Check HDF file structure."<<endl;
					exit(1);
				}

		
				//
				// Initialize the bas reader to read ccs reads as normal bas reads.
				//
		
				// Next, the location of the bases is in a non-standard group.
				ccsBasReader.baseCallsGroupName = "ConsensusBaseCalls";


				//
				// Read in the CCS fields that are the same as the base fields,
				// but in a different group.
				//

				//		ccsBasReader.OpenHDFFile(ccsBasFileName);
		
				//
				// Initialize the fields that are read.
				//
				ccsBasReader.IncludeField("Basecall");
				ccsBasReader.IncludeField("InsertionQV");
				ccsBasReader.IncludeField("DeletionQV");
				ccsBasReader.IncludeField("DeletionTag");
				ccsBasReader.IncludeField("SubstitutionQV");
				ccsBasReader.IncludeField("SubstitutionTag");
				ccsBasReader.IncludeField("QualityValue");
				//
				// Initialize this without opening a file.
				//

				ccsBasReader.Initialize(&this->rootGroup);
				//ccsBasReader.InitializeForReadingPulseInformation();
				//ccsBasReader.LoadRunInfo();
				/*
				 * Initialize pass information for reading.
				 */
				if (this->InitializeField(passesGroup, "AdapterHitAfter", adapterHitAfterArray, this->includedFields["AdapterHitAfter"]) == 0) return 0;
				if (this->InitializeField(passesGroup, "AdapterHitBefore", adapterHitBeforeArray, this->includedFields["AdapterHitBefore"]) == 0) return 0;
				if (this->InitializeField(passesGroup, "NumPasses", numPassesArray, this->includedFields["NumPasses"]) == 0) return 0;
				if (this->InitializeField(passesGroup, "PassDirection", passDirectionArray, this->includedFields["PassDirection"]) == 0) return 0;
				if (this->InitializeField(passesGroup, "PassNumBases", passNumBasesArray, this->includedFields["PassNumBases"]) == 0) return 0;
				if (this->InitializeField(passesGroup, "PassStartBase", passStartBaseArray, this->includedFields["PassStartBase"]) == 0) return 0;
				//
				// The following two fields are not critical.
				//
				this->InitializeField(passesGroup, "PassStartPulse", passStartPulseArray, this->includedFields["PassStartPulse"]);
				this->InitializeField(passesGroup, "PassNumPulses", passNumPulsesArray, this->includedFields["PassNumPulses"]);
	
				//
				// The zmw reader contains the group that hols all pass information
				//

				zmwReader.Initialize(&ccsBasReader.baseCallsGroup);
			}			
		}
		return 1;
	}

  unsigned int GetNumPasses(int readIndex) {
    unsigned int numPasses;
		numPassesArray.Read(readIndex, readIndex+1, &numPasses);
    return numPasses;
  }

	int GetNext(T_Sequence &ccsSequence) {
		//
		// Read in all ccs pass data.
		//

		ccsSequence.Free();
		int retVal = 0;
		if (this->curRead == ccsBasReader.nReads) {
			return 0;
		}
		if (this->curBasePos == ccsBasReader.nBases) {
			return 0;
		}
        try {
		numPassesArray.Read(this->curRead, this->curRead+1, &ccsSequence.numPasses);
		if (ccsSequence.numPasses > 0) {

			if (this->includedFields["AdapterHitAfter"]) {
				ccsSequence.adapterHitAfter.resize(ccsSequence.numPasses);
				adapterHitAfterArray.Read(curPassPos,  curPassPos + ccsSequence.numPasses, &ccsSequence.adapterHitAfter[0]);
			}
			if (this->includedFields["AdapterHitBefore"]) {
				ccsSequence.adapterHitBefore.resize(ccsSequence.numPasses);
				adapterHitBeforeArray.Read(curPassPos, curPassPos + ccsSequence.numPasses, &ccsSequence.adapterHitBefore[0]);
			}
			if (this->includedFields["PassDirection"]) {
				ccsSequence.passDirection.resize(ccsSequence.numPasses);
				passDirectionArray.Read(curPassPos,    curPassPos + ccsSequence.numPasses, &ccsSequence.passDirection[0]);
			}
			if (this->includedFields["PassNumBases"]) {
				ccsSequence.passNumBases.resize(ccsSequence.numPasses);
				passNumBasesArray.Read(curPassPos,     curPassPos + ccsSequence.numPasses, &ccsSequence.passNumBases[0]);
			}
			if (this->includedFields["PassStartBase"]) {
				ccsSequence.passStartBase.resize(ccsSequence.numPasses);
				passStartBaseArray.Read(curPassPos,    curPassPos + ccsSequence.numPasses, &ccsSequence.passStartBase[0]);
			}
			if (this->includedFields["PassStartPulse"]) {
				ccsSequence.passStartPulse.resize(ccsSequence.numPasses);
				passStartPulseArray.Read(curPassPos,   curPassPos + ccsSequence.numPasses, &ccsSequence.passStartPulse[0]);
			}
			if (this->includedFields["PassNumPulses"]) { 
				ccsSequence.passNumPulses.resize(ccsSequence.numPasses);
				passNumPulsesArray.Read(curPassPos,    curPassPos + ccsSequence.numPasses, &ccsSequence.passNumPulses[0]);			
			}
			curPassPos += ccsSequence.numPasses;

			// Read in the ccs bases

			retVal = ccsBasReader.GetNext((SMRTSequence&)ccsSequence);

			if (retVal == 0) {
				return 0;
			}
		}
		else {
			// advance a read in the ccs sequence without advancing positions.
			ccsBasReader.curRead++;
		}
		//
		// Regardless whether or not a ccs read was called, read the next
		// unrolled read, since an unrolled read is called for each zmw.
		//
		retVal = ((T_HDFBasReader<SMRTSequence>*)this)->GetNext(ccsSequence.unrolledRead);
        ccsSequence.zmwData = ccsSequence.unrolledRead.zmwData;
		ccsSequence.CopyTitle(ccsSequence.unrolledRead.title);
    string newTitle = string(ccsSequence.title) + string("/ccs");
    ccsSequence.CopyTitle(newTitle.c_str());
        } catch (DataSetIException e) {
            cout << "ERROR, could not read ccs data for CCS Sequence " 
                 << ccsSequence.unrolledRead.title << endl; 
            exit(1);
        }
		//		cout << "title: " << ccsSequence.title << endl;
		if (retVal == 0) {
			return 0;
		}
		else {
			return 1;
		}
	}
};


#endif
