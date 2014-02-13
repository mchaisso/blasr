#ifndef DATA_HDF_DATASET_COLLECTION_H_
#define DATA_HDF_DATASET_COLLECTION_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "HDFGroup.h"
#include "HDFData.h"


using namespace std;
using namespace std;
class DatasetCollection {
 public:
	vector<string> fieldNames;
	map<string,bool> includedFields;
	map<string,bool> requiredFields;

	void MakeFieldRequired(string &fieldName) {
		includedFields[fieldName] = true;
		requiredFields[fieldName] = true;
	}

	void MakeFieldOptional(string &fieldName) {
		includedFields[fieldName] = true;
		requiredFields[fieldName] = false;
	}

	void InitializeAllFields(bool value) {
		int f;
		for (f = 0; f < fieldNames.size(); f++ ) {
			includedFields[fieldNames[f]] = value;
		}
	}

	void InitializeFields(vector<string> &fieldList) {
		int i;
		for (i = 0; i < fieldList.size(); i++) {
			includedFields[fieldList[i]] = true;
		}
	}

	void InitializeFields(vector<char*> &fieldList) {
		int i;
		InitializeAllFields(false);
		for (i = 0; i < fieldList.size(); i++) {
			includedFields[fieldList[i]] = true;
		}
	}
	
	int IncludeField(string fieldName) {
		if (includedFields.find(fieldName) == includedFields.end()) {
			return 0;
		}	
		else {
			includedFields[fieldName] = true;
		}
        return 1;
	}

	int ExcludeField(string fieldName) {
		if (includedFields.find(fieldName) == includedFields.end()) {
			return 0;
		}	
		else {
			includedFields[fieldName] = false;
		}
		return 1;
	}
	
	bool FieldIsIncluded(string fieldName) {
		if (includedFields.find(fieldName) == includedFields.end()) {
			return false;
		}
		else {
			return includedFields[fieldName];
		}
	}
			
	bool ContainsField(string fieldName) {
		int f;
		for (f = 0; f < fieldNames.size(); f++) {
			if (fieldNames[f] == fieldName) return true;
		}
		return false;
	}

	template <typename T_Dataset>
	bool InitializeDataset(HDFGroup &group, T_Dataset &dataset, string datasetName) {
		//
		// Perform initialization of the dataset in a way that keep track
		// of which datasets in the collection are present.
		//
		if (includedFields[datasetName]) {
			if (dataset.Initialize(group, datasetName) == false) {
				if (requiredFields[datasetName]) {
					return false;
				}
				else {
					//
					// This field was supposed to be included but it either does
					// not exist or there was a problem otherwise in creating
					// it.  Don't try and read from it later on.
					//
					includedFields[datasetName] = false;
				}
			}
		}
		return true;
	}

};
	
		

#endif
