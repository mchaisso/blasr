#ifndef DATASTRUCTURES_SAF_ALN_INFO_H_
#define DATASTRUCTURES_SAF_ALN_INFO_H_

#include <stdint.h>
#include <vector>
#include "datastructures/alignment/CmpAlignment.h"
#include "Types.h"
#include <inttypes.h>

class AlnInfo {
 public:
 vector<CmpAlignment> alignments;
 UInt nAlignments;
 uint64_t lastRow;
};


#endif
