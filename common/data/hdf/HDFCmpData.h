#ifndef DATA_HDF_HDF_CMP_DATA_H_
#define DATA_HDF_HDF_CMP_DATA_H_

#include "data/hdf/HDFAtom.h"
#include "data/hdf/HDFCmpRefAlignmentGroup.h"

class HDFCmpData {
 public:
  HDFAtom<string> commandLine;
  H5File hdfCmpFile;
  HDFArray<int> movieNameIdArray;
  HDFArray<string> movieNameArray;

  HDFArray<string> readGroupPathArray;
  HDFArray<int>    readGroupPathIdArray;

  HDFArray<int>    refSeqNameIdArray;
  HDFArray<string> refSeqNameArray;
  static const int NCols=22;
  vector<HDFAtom<string> >  colNameAtoms;
  vector<HDFCmpRefAlignmentGroup*> refAlignGroups;
  map<string, int> nameToAlignmentGroupIndex;
  /*  static const char *colNameIds[] = {"00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
                              "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                              "20", "21"};*/
  static const char *colNameIds[];
  
  void Close() {
    hdfCmpFile.close();
  }
};

const char * HDFCmpData::colNameIds[] = {"00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
                                         "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                         "20", "21"};


#endif
