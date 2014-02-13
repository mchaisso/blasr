#include "data/hdf/HDFCmpFile.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFPlsReader.h"
#include "data/hdf/PlatformId.h"
#include "datastructures/alignment/CmpFile.h"
#include "datastructures/alignment/CmpAlignment.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "datastructures/reads/BaseFile.h"
#include "datastructures/reads/PulseFile.h"
#include "utils/FileOfFileNames.h"
#include <map>
#include <set>
#include <string>

typedef map<int,int> HoleToStartMap;
typedef map<string, HoleToStartMap*> MovieToHoleMap;
typedef map<string, int> MovieNameToArrayIndex;


#include "FASTAReader.h"
#include "FASTASequence.h"
#include "SMRTSequence.h"
#include <map>
#include <vector>
using namespace std;


int main(int argc, char* argv[]) {

	string cmpFileName;
	string refFileName;
	string readsFileName;
	if (argc <= 2) {
		cout << "usage: computeRefCoverage cmpfile reffile readfile" << endl;
		cout << "  computeRefCoverage Print useful stats of a cmp file" << endl;
		exit(1);
	}
	vector<int> refPositions;
	cmpFileName = argv[1];
	refFileName = argv[2];

	CmpFile cmpFile;
	FASTASequence ref;
	FASTAReader reader;

	reader.Initialize(refFileName);
	reader.GetNext(ref);

	HDFBasReader basReader;

	SMRTSequence seq, *seqPtr;

	vector<int> refCoverage;
	refCoverage.resize(ref.length);
	std::fill(refCoverage.begin(), refCoverage.end(), 0);
	/*
	 * These guys pull information from the same pls file.
	 */
	HDFCmpFile<CmpAlignment> cmpReader;


	if (cmpReader.Initialize(cmpFileName) == 0) {
		cout << "ERROR, could not open the cmp file." << endl;
		exit(1);
	}
	
	
	cmpReader.Read(cmpFile);
	UInt alignmentIndex;

	//	movieIndexSets.resize(nMovies);
	for (alignmentIndex = 0; alignmentIndex < cmpFile.alnInfo.alignments.size(); alignmentIndex++) {
        int alnGroupId  = cmpFile.alnInfo.alignments[alignmentIndex].GetAlnGroupId();
        int refGroupId  = cmpFile.alnInfo.alignments[alignmentIndex].GetRefGroupId();

        if (cmpReader.refGroupIdToArrayIndex.find(refGroupId) == cmpReader.refGroupIdToArrayIndex.end()) {
            cout << "ERROR!  An alignment " << alignmentIndex << " is specified with reference group " << endl
                 << refGroupId << " that is not found as an alignment group." << endl;
            exit(1);
        }
        int refGroupIndex = cmpReader.refGroupIdToArrayIndex[refGroupId];

        //
        // Now find the group containing the alignment.
        //
        if (cmpReader.alnGroupIdToReadGroupName.find(alnGroupId) == cmpReader.alnGroupIdToReadGroupName.end()) {
            cout << "ERROR!  An alignment " << alignmentIndex << " is specified with alignment group " << endl
                 << alnGroupId << " that is not found." << endl;
        }

        string readGroupName = cmpReader.alnGroupIdToReadGroupName[alnGroupId];
        if (cmpReader.refAlignGroups[refGroupIndex]->experimentNameToIndex.find(readGroupName) == 
            cmpReader.refAlignGroups[refGroupIndex]->experimentNameToIndex.end()) {
            cout << "ERROR!  An alignment " << alignmentIndex << " is specified with read group name " << endl
                 << readGroupName << " that is not found." << endl;
            exit(1);
        }
        int readGroupArrayIndex = cmpReader.refAlignGroups[refGroupIndex]->experimentNameToIndex[readGroupName];

		vector<char> alignedSequence, alignedTarget;

		// This read overlaps one of the ref positions.
		
		UInt offsetEnd, offsetBegin;
				
		offsetEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd();
		offsetBegin = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
		vector<unsigned char> byteAlignment;
		int alignedSequenceLength = offsetEnd - offsetBegin;
		if (alignedSequenceLength >= 0) {
			alignedSequence.resize(alignedSequenceLength);
			alignedTarget.resize(alignedSequenceLength);
			byteAlignment.resize(alignedSequenceLength);
		}

		cmpReader.refAlignGroups[refGroupIndex]->readGroups[readGroupArrayIndex]->alignmentArray.Read(offsetBegin, offsetEnd, &byteAlignment[0]);
		UInt refStart = cmpFile.alnInfo.alignments[alignmentIndex].GetRefStart();
		UInt refEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetRefEnd();
		UInt readStart= cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart();
		UInt readEnd  = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd();
		//
		// Read the alignment string.
		//
		if (refGroupIndex > 0) continue;

		//
		// Convert to something we can compare easily.
		//
		alignedSequence[alignedSequence.size()-1]= '\0';
		ByteAlignmentToQueryString(&byteAlignment[0], byteAlignment.size(), &alignedSequence[0]);
		ByteAlignmentToRefString(&byteAlignment[0], byteAlignment.size(), &alignedTarget[0]);
		int gi, i;
		gi = 0;
		int refStrand =  cmpFile.alnInfo.alignments[alignmentIndex].GetRCRefStrand();
		if (refStrand == 1) {
			// revcomp the ref strand
			vector<char> rcAlignedTarget, rcAlignedQuery;
			int t;
			rcAlignedTarget.resize(alignedTarget.size());
			rcAlignedQuery.resize(alignedSequence.size());
			for (t = 0; t < alignedTarget.size(); t++) {
				if (alignedTarget[t] == ' ') {
					rcAlignedTarget[alignedTarget.size() - t - 1] = ' ';
				}
				else {
					rcAlignedTarget[alignedTarget.size() - t - 1] = ReverseComplementNuc[alignedTarget[t]];
				}
				if (alignedSequence[t] == ' '){ 
					rcAlignedQuery[alignedTarget.size()  - t - 1] = ' ';
				}
				else {
					rcAlignedQuery[alignedTarget.size() - t - 1] = ReverseComplementNuc[alignedTarget[t]];
				}
			}
			alignedTarget = rcAlignedTarget;
			alignedSequence = rcAlignedQuery;
		}
		
		int holeNumber = cmpFile.alnInfo.alignments[alignmentIndex].GetHoleNumber();
		int ri = readStart;

		gi = refStart;

		for (i = 0; i < alignedTarget.size(); i++, gi++, ri++ ) {
			while(i < alignedTarget.size() and alignedTarget[i] == ' ') { 
				i++; 
			}
			if (alignedSequence[i] != ' ') {
				refCoverage[gi]++;
			}
		}
	} // end looping over regions

// Now compute the number of gaps.
	UInt pos;
	int numNotCovered = 0;
	for (pos = 0; pos < refCoverage.size(); pos++ ){
		if (refCoverage[pos] < 1) { numNotCovered++;}
	}
	if (numNotCovered > 100) {
		cout << "TOO Many!!!" << endl;
	}
	else {
		for (pos = 0; pos < refCoverage.size(); pos++ ){
			//		cout << refCoverage[pos] << endl;
			if (refCoverage[pos] < 1) {
				int left, right;
				left = right = -1;
				if (pos > 0) { left = refCoverage[pos-1];}
				if (pos < refCoverage.size()-1) {right = refCoverage[pos+1];}
				cout << pos << " " << left << " " << right << endl;
			}
		}
	}

}
