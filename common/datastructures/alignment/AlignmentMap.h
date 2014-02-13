#ifndef ALIGNMENT_MAP_H_
#define ALIGNMENT_MAP_H_

#include <vector>
#include <string> 

using namespace std;
class AlignmentMap {
 public:
	int qPos, tPos;
	vector<int> alignPos;
};


// Build a map of positions from (unaligned) bases to an aligned sequence
void CreateSequenceToAlignmentMap(const string & alignedSequence,
        vector<int> & baseToAlignmentMap) {
    baseToAlignmentMap.resize(alignedSequence.size());
    int alignedPos, unalignedPos;
    for (alignedPos=0, unalignedPos=0; 
         alignedPos < alignedSequence.size();
         alignedPos++) {
        if (not (alignedSequence[alignedPos] == ' ' or 
            alignedSequence[alignedPos] == '-')) {
            baseToAlignmentMap[unalignedPos] = alignedPos;
            unalignedPos++;
        }
    }
    baseToAlignmentMap.resize(unalignedPos);
}


#endif
