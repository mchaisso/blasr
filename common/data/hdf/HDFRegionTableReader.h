#ifndef DATA_HDF_HDF_REGION_TABLE_READER_H_
#define DATA_HDF_HDF_REGION_TABLE_READER_H_

#include <string>

#include "HDFFile.h"
#include "HDFArray.h"
#include "HDF2DArray.h"
#include "HDFAtom.h"
#include "PlatformId.h"

#include "../../datastructures/reads/RegionTable.h"

using namespace H5;
using namespace std;


class HDFRegionTableReader {
 public:
	HDFFile regionTableFile;
	HDFGroup pulseDataGroup;
	HDF2DArray<int> regions;
	
	HDFAtom<vector<string> > regionTypes;
	HDFAtom<vector<string> > regionDescriptions;
	HDFAtom<vector<string> > regionSources;
	HDFAtom<vector<string> > columnNames;
	int curRow;
	int nRows;
	bool fileContainsRegionTable;

	int Initialize(string &regionTableFileName, 
            const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT) {
		/*
		 * Initialize access to the HDF file.
		 */
		try {
			regionTableFile.Open(regionTableFileName.c_str(), H5F_ACC_RDONLY, fileAccPropList);
		}
		catch (Exception &e) {
			cout << e.getDetailMsg() << endl;
			return 0;
		}
		if (pulseDataGroup.Initialize(regionTableFile.rootGroup, "PulseData") == 0) {
			return 0;
		}
		if (pulseDataGroup.ContainsObject("Regions") == 0) {
			fileContainsRegionTable = false;
			return 0;
		}
		else {
			fileContainsRegionTable = true;
		}

		if (regions.Initialize(pulseDataGroup, "Regions") == 0) {
			return 0;
		}

		nRows = regions.GetNRows();

		if (columnNames.Initialize(regions.dataset, "ColumnNames") == 0) {
			return 0;
		}
		if (regionTypes.Initialize(regions.dataset, "RegionTypes") == 0) {
			return 0;
		}
		if (regionDescriptions.Initialize(regions.dataset, "RegionDescriptions") == 0) {
			return 0;
		}
		if (regionSources.Initialize(regions.dataset,  "RegionSources") == 0) {
			return 0;
		}
		
		curRow = 0;
		return 1;
	}

	int GetNext(RegionAnnotation &annotation) {
		//
		// Bail with no-op if this is the last row.
		//

		if (fileContainsRegionTable == false) {
			return 0;
		}

		if (curRow == nRows) {
			return 0;
		}

		regions.Read(curRow, curRow+1, annotation.row);
		++curRow;
		return 1;
	}	

	void RegionTypesToMap(RegionTable &table) {
		int i;
		table.regionTypeEnums.resize(table.regionTypes.size());
		for (i = 0;i < table.regionTypes.size(); i++) {
			if (table.regionTypes[i] == "GlobalAccuracy") {
				table.regionTypeEnums[i] = GlobalAccuracy;
			}
			else if (table.regionTypes[i] == "HQRegion") {
				table.regionTypeEnums[i] = HQRegion;
			}
			else if (table.regionTypes[i] == "Adapter") {
				table.regionTypeEnums[i] = Adapter;
			}
			else if (table.regionTypes[i] == "Insert") {
				table.regionTypeEnums[i] = Insert;
			}
			else if (table.regionTypes[i] == "Accuracy") {
				table.regionTypeEnums[i] = Insert;
			}
            else if (table.regionTypes[i] == "ArtifactRegion") {
				table.regionTypeEnums[i] = ArtifactRegion;
			}
			else {
				cout << "ERROR! Region Type " << table.regionTypes[i] << " is not supported.  Check Enumerations.h" << endl;
				assert(0);
			}
		}
	}

	int ReadTableAttributes(RegionTable &table) {
		if (fileContainsRegionTable == false) {
			return 0;
		}
		columnNames.Read(table.columnNames);
		regionTypes.Read(table.regionTypes);
		RegionTypesToMap(table);
		regionDescriptions.Read(table.regionDescriptions);
		regionSources.Read(table.regionSources);
    // All ok.
    return 1;
	}
	
	void Close() {
		pulseDataGroup.Close();
		regions.Close();
		regionTableFile.Close();
	}

	void ReadTable(RegionTable &table) {
		if (fileContainsRegionTable == false) {
			return;
		}
		ReadTableAttributes(table);
		table.table.resize(nRows);
		int i = 0;
		while(GetNext(table.table[curRow])) {
			i++;
		}
	}
	

};


#endif
