#ifndef DATA_HDF_HDF_REF_INFO
#define DATA_HDF_HDF_REF_INFO

#include "datastructures/saf/RefInfo.h"

class HDFRefInfoGroup {
 public:
  HDFGroup refInfoGroup;
  HDFArray<string> fullNameArray;
  HDFArray<uint32_t> idArray;
  HDFArray<uint32_t> lengthArray;
  HDFArray<string> md5Array;

  bool Create(HDFGroup &parent) {
    parent.AddGroup("RefInfo");
    if (refInfoGroup.Initialize(parent.group, "RefInfo") == 0) { return 0; }
    
    fullNameArray.Create(refInfoGroup, "FullName");
    idArray.Create(refInfoGroup, "ID");
    lengthArray.Create(refInfoGroup, "Length");
    md5Array.Create(refInfoGroup, "MD5");
    return true;
  }

  void AddRefInfo(string &fullName, unsigned int id, unsigned int length, string md5) {
    fullNameArray.Write(&fullName, 1);
    idArray.Write(&id, 1);
    lengthArray.Write(&length, 1);
    md5Array.Write(&md5, 1);
  }
  
  int Initialize(HDFGroup &parentGroup) {
    if (refInfoGroup.Initialize(parentGroup.group, "RefInfo") == 0) { return 0; }
    if (fullNameArray.Initialize(refInfoGroup, "FullName") == 0) { return 0;}
    if (idArray.Initialize(refInfoGroup,"ID") == 0) { return 0;}
    if (lengthArray.Initialize(refInfoGroup, "Length") == 0) { return 0;}
    if (md5Array.Initialize(refInfoGroup, "MD5") == 0) { return 0;}
    return 1;
  }

  ~HDFRefInfoGroup() {
    refInfoGroup.Close();
  }

  void Read(RefInfo &refInfo) {
    UInt nRow = fullNameArray.size();
    refInfo.fullName.resize(nRow);
    refInfo.id.resize(nRow);
    refInfo.length.resize(nRow);
    refInfo.md5.resize(nRow);
    /*
    if (refInfo.fullName.size() != refInfo.id.size() or
        refInfo.id.size() != refInfo.length.size() or
        refInfo.length.size() != refInfo.md5.size()) {
      cout << "Error with the RefInfo group in a cmp.h5 file.  The datasets " << endl
           << "are of different lengths but should be the same." << endl;
      exit(1);
    }
    */
    int i;
    for (i = 0; i < refInfo.fullName.size();i++) {
      fullNameArray.Read(i,i+1,&refInfo.fullName[i]);
    }
    for (i = 0; i < refInfo.md5.size(); i++) {
      md5Array.Read(i, i+1, &refInfo.md5[i]);
    }
    lengthArray.Read(0, refInfo.length.size(), &refInfo.length[0]);
    idArray.Read(0, refInfo.id.size(), &refInfo.id[0]);
  }
};


#endif
