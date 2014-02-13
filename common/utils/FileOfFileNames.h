#ifndef UTILS_FILE_OF_FILE_NAMES_H_
#define UTILS_FILE_OF_FILE_NAMES_H_
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../utils.h"
class FileOfFileNames {
 public:

	static void StoreFileOrFileList(string fileName, vector<string> &fofnList) {
		if (IsFOFN(fileName)) {
			FOFNToList(fileName, fofnList);
		}
		else {
			fofnList.push_back(fileName);
		}
	}

	static void FOFNToList(string &fofnFileName, vector<string> &fofnList) {
		ifstream fofnIn;
		CrucialOpen(fofnFileName, fofnIn);
		while(fofnIn) {
			string name;
			getline(fofnIn, name);
			if (name.size() > 0) {
				fofnList.push_back(name);
			}
		}
	}

	static bool IsFOFN(string &fileName) {
		string::size_type dotPos = fileName.rfind(".");
		bool fileNameIsFOFN = false;
		if (dotPos != string::npos) {
			string extension;
			extension.assign(fileName, dotPos+1, fileName.size() - (dotPos+1));
			if (extension == "fofn") {
					return true;
				}
		}
		return false;
	}


  static int ExpandFileNameList(vector<string> &fileNames) {
    int rfn;
    vector<string> expandedFileNames;
    for (rfn = 0; rfn < fileNames.size(); rfn++) {
      if (FileOfFileNames::IsFOFN(fileNames[rfn])) {
        string fofnFileName = fileNames[rfn];
        vector<string> fofnFileNames;
        FileOfFileNames::FOFNToList(fofnFileName, fofnFileNames);
        int f;
        for (f = 0; f < fofnFileNames.size(); f++) {
          if (FileOfFileNames::IsFOFN(fofnFileNames[f])) {
            cout << "ERROR. Nested File of File Names are not allowed. " << endl;
            exit(1);
          }
          expandedFileNames.push_back(fofnFileNames[f]);
        }
      }
      else {
        expandedFileNames.push_back(fileNames[rfn]);
      }
    }
    fileNames = expandedFileNames;
    return fileNames.size();
  }

};

#endif
