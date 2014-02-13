#ifndef DATA_HDF_HDF_ALN_INFO_GROUP_H_
#define DATA_HDF_HDF_ALN_INFO_GROUP_H_

#include "HDFGroup.h"
#include "HDFAtom.h"
#include "HDF2DArray.h"
#include <inttypes.h>

class HDFAlnInfoGroup {
 public:
  HDFGroup alnInfoGroup;
  HDF2DArray<unsigned int> alnIndexArray;
  HDFArray<float> startTime;
  static const int NCols=22;
  HDFArray<unsigned int> numPasses;
  HDFAtom<vector<string> > columnNames;


  int InitializeNumPasses() {
    numPasses.Initialize(alnInfoGroup, "NumPasses");
    return 1;
  }

  void InitializeDefaultColumnNames(vector<string> &defaultColumnNames) {
    defaultColumnNames.push_back("AlnID");
    defaultColumnNames.push_back("AlnGroupID");
    defaultColumnNames.push_back("MovieID");
    defaultColumnNames.push_back("RefGroupID");
    defaultColumnNames.push_back("tStart");
    defaultColumnNames.push_back("tEnd");
    defaultColumnNames.push_back("RCRefStrand");
    defaultColumnNames.push_back("HoleNumber");
    defaultColumnNames.push_back("SetNumber");
    defaultColumnNames.push_back("StrobeNumber");
    defaultColumnNames.push_back("MoleculeID");
    defaultColumnNames.push_back("rStart");
    defaultColumnNames.push_back("rEnd");
    defaultColumnNames.push_back("MapQV");
    defaultColumnNames.push_back("nM");
    defaultColumnNames.push_back("nMM");
    defaultColumnNames.push_back("nIns");
    defaultColumnNames.push_back("nDel");
    defaultColumnNames.push_back("Offset_begin");
    defaultColumnNames.push_back("Offset_end");
    defaultColumnNames.push_back("nBackRead");
    defaultColumnNames.push_back("nReadOverlap");
  }
  
  bool Create(HDFGroup &parent) {
    parent.AddGroup("AlnInfo");
    // Make sure it was created, and intialize this group to reference the newly created one.
    if (alnInfoGroup.Initialize(parent.group, "AlnInfo") == 0) { return 0; }
    vector<string> defaultColumnNames;
    InitializeDefaultColumnNames(defaultColumnNames);
    columnNames.Create(alnInfoGroup.group, "ColumnNames", defaultColumnNames);
    
    alnIndexArray.Create(&alnInfoGroup.group, "AlnIndex", defaultColumnNames.size());
    return true;
  }
  int Initialize(HDFGroup &rootGroup) {
    if (alnInfoGroup.Initialize(rootGroup.group, "AlnInfo") == 0) { return 0; }
    if (alnIndexArray.Initialize(alnInfoGroup, "AlnIndex") == 0) { return 0; }
    /*
     * This functionality should go into the python.
    if (!alnIndexArray.ContainsAttribute("ColumnNames")) {
      try {
        vector<string> defaultColumnNames;
        InitializeDefaultColumnNames(defaultColumnNames);
        columnNames.Create(alnIndexArray.dataset, "ColumnNames", defaultColumnNames);
      }
      catch(Execption e) {
        //
        // If the dataset is not writable
      }
      }
    */
    return 1;
  }

  HDFAtom<int> frameRate;
  
  ~HDFAlnInfoGroup() {
    alnInfoGroup.Close();
  }
  
  // Return size of /AlnInfo/AlnIndex in KB
  UInt GetAlnIndexSize() {
    return alnIndexArray.GetNRows() / 1024 * sizeof (unsigned int) * NCols;
  }

  void Read(AlnInfo &alnInfo) {

    int nAlignments = alnIndexArray.GetNRows();
    alnInfo.alignments.resize(nAlignments);
    UInt alignmentIndex;
    UInt alignmentRow[NCols];
    for (alignmentIndex = 0; alignmentIndex < nAlignments; alignmentIndex++) {
      // Input the values.
      alnIndexArray.Read(alignmentIndex, alignmentIndex + 1, alignmentRow);
      alnInfo.alignments[alignmentIndex].StoreAlignmentIndex(alignmentRow, NCols);
    }
  }

  int GetNAlignments() {
    return alnIndexArray.GetNRows();
  }

  unsigned int WriteAlnIndex(vector<unsigned int> &aln) {
    alnIndexArray.WriteRow(&aln[0], aln.size());
    return alnIndexArray.GetNRows();
  }

  void ReadCmpAlignment(UInt alignmentIndex, CmpAlignment &cmpAlignment) {
    UInt alignmentRow[NCols];
    alnIndexArray.Read(alignmentIndex, alignmentIndex + 1, alignmentRow);
    cmpAlignment.StoreAlignmentIndex(alignmentRow, NCols);
  }
};

#endif
