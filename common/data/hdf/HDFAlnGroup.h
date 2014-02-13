#ifndef DATA_HDF_HDF_ALN_GROUP_H_
#define DATA_HDF_HDF_ALN_GROUP_H_

#include "HDFArray.h"
#include "HDFGroup.h"

class HDFAlnGroup {
 public:
  HDFGroup alnGroup;
  HDFArray<unsigned int> idArray;
  HDFArray<string> pathArray;

  void Initialize(HDFGroup &parent) {
    alnGroup.Initialize(parent.group, "AlnGroup");
    idArray.Initialize(alnGroup.group, "ID");
  }
  
  void Read(AlnGroup &aln) {
    // Seem to write data in this HDFAlnGroup obj to AlnGroup & aln
    int idNElem = idArray.arrayLength;
    int pathNElem = pathArray.arrayLength;
    if (idNElem > 0) {

      aln.id.resize(idNElem);
      idArray.Read(0, idNElem);

      aln.path.resize(pathNElem);
      unsigned int i;
      for (i = 0; i < pathNElem; i++) {
        pathArray.Read(i, i+1, &aln.path[i]);
      }
    }
  }

  int AddPath(string &path) {
    int id;
    pathArray.Write(&path, 1);
    id = pathArray.size();
    idArray.Write(&id, 1);
    return id;
  }

  ~HDFAlnGroup() {
    alnGroup.Close();
  }

};


#endif
