#ifndef DATA_HDF_HDF_SMRT_SEQUENCE_READER
#define DATA_HDF_HDF_SMRT_SEQUENCE_READER

#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFZMWReader.h"
#include "SMRTSequence.h"

template<typename T_SMRT_Sequence>
class HDFSMRTSequenceReader : public HDFBasReader {
 public:
	HDFZMWReader zmwReader;
	bool readQuality;
	int Initialize(string hdfBasFileName, bool _readQuality=true, const H5::FileAccPropList & fileAccProplist=H5::FileAccPropList::DEFAULT) {
		HDFBasReader::Initialize(hdfBasFileName, fileAccPropList);
		zmwReader.Initialize(hdfBasFile);
		readQuality = _readQuality;
		if (baseCallsGroup.ContainsObject("WidthInFrames")) {
			useWidthInFrames = InitializeField(baseCallsGroup, "WidthInFrames");
		}
		if (baseCallsGroup.ContainsObject("PreBaseFrames")) {
			usePreBaseFrames = InitializeField(baseCallsGroup, "PreBaseFrames");
		}
		if (baseCallsGroup.ContainsObject("PulseIndex")) {
			usePulseIndex = InitializeField(baseCallsGroup, "PulseIndex");
		}
		
	}


	int GetNext(T_SMRT_Sequence &seq) {
		int retVal;
		//
		// Copy the curBasePos from the bas reader since it gets advanced
		// in the GetNext function.
		//
		DNALength curBasePosCopy = curBasePos;
		retVal = HDFBasReader::GetNext(seq);
		//
		// Bail now if the file is already done
		if (retVal == 0) {
			return retVal;
		}
		zmwReader.GetNext(seq.zmwData);
		return retVal;
	}
	int Advance(int nSteps) {
		int retVal;
		retVal = HDFBasReader::Advance(nSteps);
		zmwReader.Advance(nSteps);
		return retVal;
	}


};

template<>
int HDFSMRTSequenceReader<FASTASequence>::GetNext(FASTASequence &seq) {
		int retVal;
		if (readQuality) {
			retVal = HDFBasReader::GetNext(seq);
		}
		//
		// Bail now if the file is already done
		if (retVal == 0) {
			return retVal;
		}
		zmwReader.GetNext(seq.zmwData);
		return retVal;
	}

#endif
