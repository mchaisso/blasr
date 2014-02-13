#ifndef DATA_HDF_SENTINAL_FILE
#define DATA_HDF_SENTINAL_FILE

#include "HDFArray.h"
#include "HDF2DArray.h"

class HDFSentinal {
 public:
  HDFArray<string> parts;
  HDF2DArray<unsigned int> holeLookup;
  HDFGroup multiPartGroup;

  void Initialize(HDFGroup &rootGroup) {
    multiPartGroup.Initialize(rootGroup, "MultiPart");
    parts.Initialize(multiPartGroup, "Parts");
    holeLookup.Initialize(multiPartGroup, "HoleLookup");
  }
};


#endif
