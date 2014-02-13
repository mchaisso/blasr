#ifndef DATA_HDF_HDF_DATA_H_
#define DATA_HDF_HDF_DATA_H_


#include <assert.h>
#include <iostream>
#include <vector>
#include <string>
#include "H5Cpp.h"
#include "HDFGroup.h"
#include "HDFConfig.h"
#include "HDFAttributable.h"

using namespace std;
using namespace H5;

class HDFData : public HDFAttributable {
 public:
	DataSet   dataset;
	DataSpace dataspace;
	DataSpace sourceSpace, destSpace;
	DataSpace fullSourceSpace;
	bool      fileDataSpaceInitialized;
	CommonFG  *container;
	string    datasetName;
	bool      isInitialized;


  H5Object* GetObject() {
    return &dataset;
  }
	
  HDFData(CommonFG* _container, string _datasetName) {
		container   = _container;
		datasetName = _datasetName;
		fileDataSpaceInitialized = false;
		isInitialized = false;
	}

	HDFData() {
    container = NULL;
		fileDataSpaceInitialized = false;
		isInitialized = false;
	}

	bool IsInitialized() {
		return isInitialized;
	}

  //
  // Allow derived classes to be initialized generically.
  //
  virtual int Initialize(HDFGroup &parentGroup, const string &datasetName) { 
    cout << "ERROR! Only a subclass should call this." << endl;
    exit(1);
  }

  int BaseInitializeDataset(CommonFG &hdfFile, string _datasetName) {
    dataset   = hdfFile.openDataSet(_datasetName.c_str());
    isInitialized = true;
    fileDataSpaceInitialized = true;
    return 1;
  }
    
  int InitializeDataset(HDFGroup &group, string _datasetName) {
    return InitializeDataset(group.group, _datasetName);
  }

	int InitializeDataset(CommonFG &hdfFile, string _datasetName) {
		try {
			datasetName = _datasetName;
			dataset   = hdfFile.openDataSet(_datasetName.c_str());
			isInitialized = true;
			fileDataSpaceInitialized = true;
		}
    catch(FileIException &e) {
			cerr << e.getDetailMsg() <<endl;
			return 0;
		}
		catch(GroupIException &e) {
			cerr << e.getDetailMsg() << endl;
			return 0;
		}
		catch(H5::Exception &e) {
			cerr << e.getDetailMsg() << endl;
			return 0;
		}
		return 1;
	}

	void Close() {
        if (isInitialized) {
            dataspace.close();
		    dataset.close();
            isInitialized = false;
        }
	}
};

#endif
