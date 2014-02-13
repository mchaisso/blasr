#ifndef DATASTRUCTURES_ALIGNMENT_CMP_FILE_H_
#define DATASTRUCTURES_ALIGNMENT_CMP_FILE_H_

#include "datastructures/reads/ReadType.h"
#include "datastructures/alignment/CmpIndexedStringTable.h"
#include "datastructures/saf/AlnGroup.h"
#include "datastructures/saf/AlnInfo.h"
#include "datastructures/saf/RefGroup.h"
#include "datastructures/saf/RefInfo.h"
#include "datastructures/saf/MovieInfo.h"
#include "Enumerations.h"

#include <vector>
using namespace std;

class CmpFile {
 public:
	int lastRow;
	string readTypeString, index, version, commandLine;
  ReadType::ReadTypeEnum readType;

  void StoreReadType(string &readTypeStringP) {
    readTypeString = readTypeStringP;
    readType = ParseReadType(readTypeString);
  }

	CmpIndexedStringTable readGroupTable, movieNameTable, refSeqTable;
	vector<string> colNames;
	PlatformType platformId;
	AlnGroup  alnGroup;
	AlnInfo   alnInfo;
	RefGroup  refGroup;
	RefInfo   refInfo;
	MovieInfo movieInfo;
};

#endif
