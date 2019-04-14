#ifndef ALIGNMENT_READER_AGGLOMERATE
#define ALIGNMENT_READER_AGGLOMERATE

#include <cstdlib>

#include "BaseSequenceIO.h"

#include "FASTAReader.h"
#include "FASTQReader.h"
#include "CCSSequence.h"
#include "SMRTSequence.h"
#include "Enumerations.h"
//#include "data/hdf/HDFBasReader.h"
//#include "data/hdf/HDFCCSReader.h"
#include "utils/StringUtils.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "datastructures/alignmentset/SAMToSMRTSequence.h"
#include "datastructures/reads/ScanData.h"


class ReaderAgglomerate : public BaseSequenceIO {
	FASTAReader fastaReader;
	FASTQReader fastqReader;
	
	int readQuality;
	int stride;
	int start;
	float subsample;
	bool useRegionTable;
	bool ignoreCCS;
	int lineNumber;
	bool initialized;
	int holeNumberI;

 public:
	vector<int>    holeNumbers;
	bool IsInitialized() {
		return initialized;
	}
	bool UseHoleNumbers() {
		return holeNumbers.size() > 0;
	}

	//
	// Create containers for reading hdf
	//
	/*
	T_HDFBasReader<SMRTSequence>  hdfBasReader;
	HDFCCSReader<CCSSequence>     hdfCcsReader;
	*/
	vector<SMRTSequence>          readBuffer;
	vector<CCSSequence>           ccsBuffer;


	//
	// Container for reading sam
	//
  SAMReader<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> samReader;
	
  string readGroupId;


  void SetToUpper() {
    fastaReader.SetToUpper();
  }
	void InitializeParameters() {
		start  = 0;
		stride = 1;
		subsample = 1.1;
		readQuality = 1;
		useRegionTable = true;
		ignoreCCS = true;
		initialized = false;
	}

	ReaderAgglomerate() {
		
		InitializeParameters();
	}
	
	ReaderAgglomerate(float _subsample) {
		this->InitializeParameters();
		subsample = _subsample;
	}
	 
	ReaderAgglomerate(int _stride) {
		this->InitializeParameters();
		stride = _stride;
	}
	
	ReaderAgglomerate(int _start, int _stride) {
		this->InitializeParameters();
		start  = _start;
		stride = _stride;
	}

	void GetScanData(ScanData &scanData) {
		if (fileType == Fasta || fileType == Fastq) {
		  scanData.Clear();
		  scanData.movieName = fileName;
		}
		else if (fileType == HDFPulse || fileType == HDFBase || fileType == HDFCCS) {
		  cerr << "HDF input is no longer supported" << endl;
		  exit(1);
		}
	}

	void GetChangelistId(string &softwareVersion) {
		if (fileType == Fasta || fileType == Fastq) {
			softwareVersion = "";
		}
		else if (fileType == HDFPulse || fileType == HDFBase || fileType == HDFCCS) {
		  cerr << "HDF input is no longer supported" << endl;
		  exit(1);
		  //			hdfBasReader.GetChangelistId(softwareVersion);
		}
	}

  void GetMovieName(string &movieName) {
    if (fileType == Fasta || fileType == Fastq) {
      movieName = fileName;
    }
    else if (fileType == HDFPulse || fileType == HDFBase || fileType == HDFCCS) {
      cerr << "HDF input is no longer supported" << endl;
      exit(1);
      //            movieName = hdfBasReader.GetMovieName();
    }
  }

	bool FileHasZMWInformation() {
	  return false;
	  //		return (fileType == HDFPulse || fileType == HDFBase || fileType == HDFCCS);
	}

	void SkipReadQuality() {
		readQuality = 0;
	}
	
	void IgnoreCCS() {
		ignoreCCS = true;
	}
	
	void UseCCS() {
		ignoreCCS = false;
		//hdfBasReader.SetReadBasesFromCCS();
	}

	int Initialize(string &pFileName, string fileTypeSuffix="") {
		if ((fileTypeSuffix == "" and DetermineFileTypeByExtension(pFileName, fileType, false)) or
				DetermineFileTypeByExtension("." + fileTypeSuffix, fileType, false)) {
			fileName = pFileName;

			return Initialize();
		}
		else {
			return false;
		}
	}
	
	bool SetReadFileName(string &pFileName, string fileTypeSuffix="") {
		if ((fileTypeSuffix == "" and DetermineFileTypeByExtension(pFileName, fileType, false)) or
				DetermineFileTypeByExtension("." + fileTypeSuffix, fileType, false)) {
			fileName = pFileName;
			return true;
		}
		else {
			return false;
		}
	}
		
	int Initialize(FileType &pFileType, string &pFileName) {
		SetFiles(pFileType, pFileName);
		
		return Initialize();
	}

	bool HasRegionTable() {
	  return false;
	}

	void InitializeHoleNumbers(vector<int> &hn) {
		holeNumbers = hn;
		holeNumberI = 0;
	}

	int Initialize() {
		int init = 1;
		switch(fileType) {
		case Fasta:
			init = fastaReader.Init(fileName);
			break;
		case Fastq:
			init = fastqReader.Init(fileName);
			break;
		case HDFPulse:
		case HDFBase:
		  break;
		case SAM_READ:
			init = samReader.Initialize(fileName);
			AlignmentSet<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> alignmentSet;
			samReader.ReadHeader(alignmentSet);
			break;
		}
    readGroupId = "";
		if (init == 0 || (start > 0 && Advance(start) == 0) ){
			return 0;
		};

    string movieName;
    GetMovieName(movieName);
    MakeMD5(movieName, readGroupId, 10);
		initialized = true;    
		return 1;
	}

	ReaderAgglomerate &operator=(ReaderAgglomerate &rhs) {
		fileType     = rhs.fileType;
		fileName = rhs.fileName;
		return *this;
	}
	
	bool Subsample(float rate) {
		bool retVal = true;
		while( (rand() % 100 + 1) > (rate * 100) and (retVal = Advance(1)));
		return retVal;
	}

	int GetNext(FASTASequence &seq) {
		int numRecords = 0;
		if (Subsample(subsample) == 0) {
			return 0;
		}
		switch(fileType) {
		case Fasta:
			numRecords = fastaReader.GetNext(seq);
			break;
		case Fastq:
			cout <<" Getting fatstq " << endl;
			numRecords = fastqReader.GetNext(seq);
			break;
		case HDFPulse:
		  cerr << "HDF no longer supported" << endl;
		  exit(1);
			break;
		}
    seq.CleanupOnFree();
		return numRecords;
	}

	int GetNext(FASTQSequence &seq) {
		int numRecords = 0;
		if (Subsample(subsample) == 0) {
			return 0;
		}
		switch(fileType) {
		case Fasta:
			numRecords = fastaReader.GetNext(seq);
			break;
		case Fastq:
			numRecords = fastqReader.GetNext(seq);
			break;
		case HDFPulse:
                case HDFBase:
		  cerr << "HDF no longer supported" << endl;
		  exit(1);

			break;
		case HDFCCS:
		  cerr << "HDF no longer supported" << endl;
		  exit(1);

			break;
		}
		if (stride > 1)
			Advance(stride-1);
		return numRecords;
	}

	int GetNext(SMRTSequence &seq) {
		int numRecords = 0;

		if (UseHoleNumbers()) {
			if (holeNumberI >= holeNumbers.size()){ 
				return 0;
			}
		}

		if (Subsample(subsample) == 0) {
			return 0;
		}
		switch(fileType) {
		case Fasta:
			numRecords = fastaReader.GetNext(seq);
			break;
		case Fastq:
			numRecords = fastqReader.GetNext(seq);
			break;
		case HDFPulse:
		case HDFBase:
		  cerr << "HDF no longer supported" << endl;
		  exit(1);
		  break;
    case HDFCCS:
      		  cerr << "HDF no longer supported" << endl;
		  exit(1);
			break;

		case SAM_READ:
			SAMAlignment samAlignment;
			lineNumber = 0;
			//
			// Read the SAM header if it exists, but it is not used.
			//
			if (samReader.GetNextAlignment(samAlignment, true) == false) {
				return 0;
			}
			ConvertSAMToSMRTSequence(samAlignment, seq);
			numRecords = 1;
			break;
		}
		if (stride > 1)
			Advance(stride-1);
		return numRecords;
	}

	int GetNext(CCSSequence &seq) {
		int numRecords = 0;
		if (Subsample(subsample) == 0) {
			return 0;
		}

		switch(fileType) {
		case Fasta:
			// This just reads in the fasta sequence as if it were a ccs sequence
			numRecords = fastaReader.GetNext(seq);
			seq.subreadStart = 0;
			seq.subreadEnd   = 0;
			break;
		case Fastq:
			numRecords = fastqReader.GetNext(seq);
			seq.subreadStart = 0;
			seq.subreadEnd   = 0;
			break;
		case HDFPulse:
    case HDFBase:
      		  cerr << "HDF no longer supported" << endl;
		  exit(1);


			break;
    case HDFCCS:
		  cerr << "HDF no longer supported" << endl;
		  exit(1);

			break;
		}

		if (stride > 1)
			Advance(stride-1);
		return numRecords;
	}

	int Advance(int nSteps) {
    int i;
		switch(fileType) {
		case Fasta:
			return fastaReader.Advance(nSteps);
		case HDFPulse:
		case HDFBase:
		case HDFCCS:
		  cerr << "HDF no longer supported" << endl;
		  exit(1);
		  break;
		case Fastq:
			return fastqReader.Advance(nSteps);
		}
    return false;
	}
	
	void Close() {
		switch(fileType) {

		case Fasta:
		case Fastq:
			fastaReader.Close();
			break;
		case HDFPulse:
		case HDFBase:
		case HDFCCS:
		  cerr << "HDF no longer supported" << endl;
		  exit(1);		  
			break;
		case SAM_READ:
			samReader.Close();
			break;
		}
		initialized = false;
	}
};

template<typename T_Sequence>
int ReadChunkByNReads(ReaderAgglomerate &reader, vector<T_Sequence> &reads, int maxNReads) {
	T_Sequence seq;
	int nReads = 0;
	while(nReads < maxNReads) {
		if (reader.GetNext(seq)) {
			reads.push_back(seq);
			++nReads;
		}
		else {
			break;
		}
	}
	return nReads;
}

template<typename T_Sequence>
int ReadChunkBySize (ReaderAgglomerate &reader, vector<T_Sequence> &reads, int maxMemorySize) {
	T_Sequence seq;
	int nReads = 0;
	int totalStorage = 0;
	while (totalStorage < maxMemorySize) {
		if (reader.GetNext(seq)) {
			reads.push_back(seq);
			totalStorage += seq.GetStorageSize();
			nReads++;
		}
		else {
			break;
		}
	}
	return nReads;
}
	 
		
			 
#endif
