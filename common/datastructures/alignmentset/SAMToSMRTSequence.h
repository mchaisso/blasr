#ifndef SAM_TO_SMRT_SEQUENCE_H_
#define SAM_TO_SMRT_SEQUENCE_H_

#include "datastructures/alignmentset/SAMAlignment.h"
#include "SMRTSequence.h"
#include <iostream>
#include <string>


int CopyFieldValues(unsigned char* &dest, string src) {
	if (src.size() == 0) {
		if (dest !=  NULL) {
			delete [] dest;
			dest = NULL;
		}
		return 0;
	}
	else {
		memcpy(dest, src.c_str(), src.size());
		return src.size();
	}
}



int AssignData( unsigned char* &field, string src) {
	if (field == NULL) {
		if (src.size() > 0) {
			cout << "ERROR. Copying into unallocated field." << endl;
			assert(0);
		}
		else {
			return 0;
		}
	}
	int res = CopyFieldValues(field, src);
	return res;
}

int AssignQualityData( unsigned char* &field, string src) {
	int res = AssignData(field, src);

	if (res == 0) {
		field = NULL;
	}		
	else {
		QualityStringToStored(field, src.size());
	}
	return res;
}

bool ConvertSAMToSMRTSequence(SAMAlignment &samAlignment, SMRTSequence &read) {
	//
	// Populate sequence / quality information.  Any empty qualities
	// are skipped.
	//
	
	read.Allocate(samAlignment.seq.size());
	read.length = samAlignment.seq.size();
	AssignData(read.seq, samAlignment.seq);
	if (samAlignment.qual != "*") {
		AssignQualityData(read.qual.data, samAlignment.qual);
	}
	else {
		read.qual.Free();
	}

	AssignQualityData(read.deletionQV.data, samAlignment.qd);
	AssignQualityData(read.insertionQV.data, samAlignment.qi);
	AssignQualityData(read.substitutionQV.data, samAlignment.qs);
	AssignQualityData(read.mergeQV.data, samAlignment.qm);
	AssignData((unsigned char*&) read.substitutionTag, samAlignment.ts);
	AssignData((unsigned char*&) read.deletionTag, samAlignment.td);

	
	read.StoreHoleNumber(samAlignment.zmw);
	read.StoreHoleStatus(0);
	read.zmwData.numEvents = read.length;
	
	memcpy(read.snr, samAlignment.snr, sizeof(float)*4);
	read.numPasses = samAlignment.numPasses;
	read.holeNumber = samAlignment.zmw;
	read.subreadStart = 0;
	read.subreadEnd   = read.length;
	read.qs = samAlignment.qStart;
	read.qe = samAlignment.qEnd;
	read.readQuality = samAlignment.readQuality;
	read.CopyTitle(samAlignment.qName);

}
#endif
