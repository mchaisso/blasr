#ifndef DATA_HDF_HDF_CMP_ROOT_GROUP_H_
#define DATA_HDF_HDF_CMP_ROOT_GROUP_H_

#include "HDFAtom.h"
#include "HDF2DArray.h"
#include "datastructures/alignment/CmpFile.h"
template <typename T_Alignment>
class HDFCmpRootGroup {
 public:
	HDFGroup rootGroup;
	HDFAtom<string> version;
	HDFAtom<string> index;
	HDFAtom<string> readType;
	HDFAtom<string> commandLine;
	HDF2DArray<string> fileLog;

	~HDFCmpRootGroup() {
		rootGroup.Close();
	}

	int Initialize(H5File &cmpFile) {
		if (rootGroup.Initialize(cmpFile, "/") == 0) { return 0; }
		if (rootGroup.ContainsObject("Version")) {
			if (version.Initialize(rootGroup.group, "Version") == 0) { return 0; }		
		}
		if (rootGroup.ContainsObject("Index")) {
			if (index.Initialize(rootGroup.group, "Index") == 0) { return 0; }
		}
		if (rootGroup.ContainsObject("ReadType")) {
			if (readType.Initialize(rootGroup.group, "ReadType") == 0) { return 0; }
		}
		if (rootGroup.ContainsObject("CommandLine")) {
			if (commandLine.Initialize(rootGroup.group, "CommandLine") == 0) { return 0; }
		}
		
		//
		// For now, disable file log initialization until
		// hdf2darray<string> is tested more thoroughly 
		//
		// if (fileLog.Initialize(rootGroup.group, "FileLog") == 0) {
		// return 0;}
		//
		return 1;
	}
	
	void ReadAttributes(CmpFile &cmpFile) {
		if (rootGroup.ContainsObject("Version")) {
		version.Read(cmpFile.version);
		}
		if (rootGroup.ContainsObject("Index")) {
			index.Read(cmpFile.index);
		}
		if (rootGroup.ContainsObject("ReadType")) {
      string readTypeString;
      readType.Read(readTypeString);
      cmpFile.StoreReadType(readTypeString);
		}
		if (rootGroup.ContainsObject("CommandLine")) {
		commandLine.Read(cmpFile.commandLine);
		}
	}
};
 
	

#endif
