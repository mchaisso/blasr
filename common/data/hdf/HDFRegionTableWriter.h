#ifndef DATA_HDF_HDF_REGION_TABLE_WRITER_H_
#define DATA_HDF_HDF_REGION_TABLE_WRITER_H_

#include <string>

#include "HDFFile.h"
#include "HDFArray.h"
#include "HDF2DArray.h"
#include "HDFAtom.h"
#include "PlatformId.h"

#include "datastructures/reads/RegionTable.h"


using namespace H5;
using namespace std;


class HDFRegionTableWriter {
 public:
	HDFFile regionTableFile;

	HDF2DArray<int> regions;
	HDFAtom<vector<string> > regionTypes;
	HDFAtom<vector<string> > regionDescriptions;
	HDFAtom<vector<string> > regionSources;
	HDFAtom<vector<string> > columnNames;
	int curRow;
	int nRows;
	HDFGroup *parentGroupPtr, pulseDataGroup;

  void CreateGroupStructure() {

    regionTableFile.rootGroup.AddGroup("PulseData");
    if (pulseDataGroup.Initialize(regionTableFile.rootGroup, "PulseData") == 0) {
      cout << "Could not create group PulseData. This is a bug." << endl;
      exit(1);
    }
    parentGroupPtr = &pulseDataGroup;
    regions.Initialize(pulseDataGroup, "Regions", RegionAnnotation::NCOLS);

  }

  int Create(string fileName) {
    H5File newFile(fileName.c_str(), H5F_ACC_TRUNC, FileCreatPropList::DEFAULT, FileAccPropList::DEFAULT);  
    regionTableFile.hdfFile.openFile(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);    
    regionTableFile.rootGroup.Initialize(regionTableFile.hdfFile, "/");
    CreateGroupStructure();
  }
  
  HDFRegionTableWriter() {
    curRow = 0;
    nRows  = 0;
    parentGroupPtr = NULL;
  }

  int Initialize(HDFGroup &parentGroupP) {
    /*
     * Initialize in an existing file.
     */
    parentGroupPtr = &parentGroupP;
    return Initialize();
  }

	int Initialize(string &regionTableFileName, 
                   const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT) {
		/*
		 * Initialize access to the HDF file.
		 */
		try {
			regionTableFile.Open(regionTableFileName, H5F_ACC_TRUNC, fileAccPropList);
		}
		catch (Exception &e) {
			cout << e.getDetailMsg() << endl;
			return 0;
		}
    CreateGroupStructure();
    return 1;
  }

  int Initialize() {
    
		if (regions.Initialize(*parentGroupPtr, "Regions", RegionAnnotation::NCOLS) == 0) {
			return 0;
		}
		nRows = regions.GetNRows();


		curRow = 0;
		return 1;
	}

  void Finalize(vector<string> &columnNamesVect, 
                 vector<string> &regionTypesVect,
                 vector<string> &regionDescriptionsVect,
                 vector<string> &regionSourcesVect) {
    //
    // Make sure data has been written to the dataset.  If not, the
    // dataset will not exist and creating the attributes will fail.
    //
    if (curRow > 0) {
      regionTypes.Create(regions.dataset, "RegionTypes", regionTypesVect);
      columnNames.Create(regions.dataset, "ColumnNames", columnNamesVect);
      regionDescriptions.Create(regions.dataset, "RegionDescriptions", regionDescriptionsVect);
      regionSources.Create(regions.dataset,  "RegionSources", regionSourcesVect);
    }
  }
  void WriteRows(vector<RegionAnnotation> &annotations) {
    int i;
    for (i = 0; i < annotations.size(); i++) {
      Write(annotations[i]);
    }
  }

	int Write(RegionAnnotation &annotation, int at=-1) {
    if (at == -1) {
      regions.WriteRow(annotation.row, annotation.NCOLS);
      ++curRow;
    }
    else {
      regions.WriteRow(annotation.row, annotation.NCOLS, at);
    }
		return 1;
	}	

  void Close() {
    regions.Close();
    //
    // Check to see if this is external
    //
    if (parentGroupPtr == &pulseDataGroup) {
      pulseDataGroup.Close();
    }
    regionTableFile.rootGroup.Close();
  }

};


#endif
