#ifndef DATA_HDF_HDF_FILE_H_
#define DATA_HDF_HDF_FILE_H_

#include <iostream>
#include <string>
#include <vector>
#include "H5Cpp.h"

#include "HDFConfig.h"
#include "HDFGroup.h"

using namespace H5;
using namespace std;


class HDFFile {
 public:
	// Make this public for easy access.
	H5File hdfFile;	
  HDFGroup rootGroup;

	HDFFile() {
	}

  //
  //  Open a file.  By default, if the file already exists, open it in
  //  read/write mode.  The only other flag that is allowed is
  //  H5F_ACC_TRUNC, which will truncate the file to zero size.
  //
  void Open(string fileName, unsigned int flags=H5F_ACC_RDWR, const FileAccPropList & fileAccPropList=FileAccPropList::DEFAULT) {
    assert (flags == H5F_ACC_RDWR || flags == H5F_ACC_TRUNC || flags == H5F_ACC_RDONLY);
    ifstream testIn(fileName.c_str());
    bool fileExists = testIn;
    bool flagsIsNotTrunc = flags != H5F_ACC_TRUNC;

    if (fileExists and H5File::isHdf5(fileName.c_str()) and flagsIsNotTrunc) {
      try {
        hdfFile.openFile(fileName.c_str(), flags, fileAccPropList);
      }
      catch (FileIException e) {
        cout << "Error opening file " << fileName << endl;
        exit(1);
      }
    }
    else {
      try {
        //
        // Open a new file with TRUNC permissions, always read/write.
        //
        FileCreatPropList filePropList;
        hsize_t ub = filePropList.getUserblock();
        filePropList.setUserblock(512);
        hdfFile = H5File(fileName.c_str(), H5F_ACC_TRUNC, filePropList);
      }
      catch (FileIException fileException) {
        cout << "Error creating file " << fileName << endl;
        exit(1);
      }
    }
    if (rootGroup.Initialize(hdfFile, "/") != 1) {
        cout << "Error initializing the root group for file " << fileName << endl;
        exit(1);
    }
  }

	void Close() {
    hdfFile.close();
	}

};


#endif
