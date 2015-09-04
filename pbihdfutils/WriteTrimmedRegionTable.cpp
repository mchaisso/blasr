#include "CommandLineParser.h"
#include "data/hdf/HDFPlsReader.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "datastructures/reads/RegionTable.h"
#include "datastructures/reads/ReadInterval.h"
#include "data/hdf/HDFBasWriter.h"
#include "utils/StringUtils.h"
#include "Enumerations.h"
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]) {


	string baxFileName, adapterAlignmentFile;
  CommandLineParser clp;
  clp.SetProgramName("writeTrimmedRegionTable");
  clp.RegisterStringOption("file.bax.h5", &baxFileName, "Input bax.h5 file.", true);
  clp.RegisterStringOption("alignments.m1", &adapterAlignmentFile, "Input bax.h5 file.", true);
	//  clp.RegisterStringOption("output.rgn.h5", &rgnFileName, "Output region file.", true);	
  clp.RegisterPreviousFlagsAsHidden();
  clp.ParseCommandLine(argc, argv);


	HDFRegionTableReader hdfRegionReader;
	HDFBasReader basReader;
	basReader.Initialize(baxFileName);
	hdfRegionReader.Initialize(baxFileName);

	RegionTable regionTable;
	hdfRegionReader.ReadTable(regionTable);
	regionTable.SortTableByHoleNumber();

	
	ifstream m1File;
	CrucialOpen(adapterAlignmentFile, m1File, std::ios::in);
  HDFRegionTableWriter regionWriter;
	HDFBasWriter writer;
	//	writer.Initialize(baxFileName, basReader.GetMovieName(), basReader.GetRunCode());
	regionWriter.Initialize(basReader.pulseDataGroup);
	while (m1File) {
		string readName, targetName;
		int qStrand, tStrand, alnScore, tStart, tEnd, tLength, qStart, qEnd, qLength, tmp;
		float ident;
		if ((m1File >> readName >> targetName >> tStrand >> qStrand >> alnScore >> ident >> tStart >> tEnd >> tLength >> qStart >> qEnd >> qLength >> tmp) == 0) {
			break;
		}
		
		vector<string> readNameParts;
		ParseSeparatedList(readName, readNameParts, '/');
		vector<string> subread;
		ParseSeparatedList(readNameParts[2], subread, '_');
		int subreadStart,  subreadEnd;
		subreadStart = atoi(subread[0].c_str());
		subreadEnd   = atoi(subread[1].c_str());
		
		int holeNumber = atoi(readNameParts[1].c_str());

		int insertStart, insertEnd;
		if (subreadEnd - qEnd > qStart - subreadStart) {
			insertStart = qEnd;
			insertEnd  = subreadEnd;
		}
		else {
			insertStart = subreadStart;
			insertEnd   = qStart;
		}
			
		int low, high;
		regionTable.LookupRegionsByHoleNumber(holeNumber, low, high);
		int i;
		for (i = low; i < high; i++) {
			RegionType regionType = regionTable.GetType(i);
			if (regionType == Insert and (regionTable.GetStart(i) == subreadStart or regionTable.GetEnd(i) == subreadEnd)) {
				regionTable.SetStart(i, insertStart);
				regionTable.SetEnd(i, insertEnd);
			}
		}
	}
	int regionIndex;
	for (regionIndex = 0; regionIndex < regionTable.table.size(); regionIndex++) {
		regionWriter.Write(regionTable.table[regionIndex], regionIndex);
	}
	/*	regionWriter.Finalize(regionTable.columnNames,
												regionTable.regionTypes, 
												regionTable.regionDescriptions, 
												regionTable.regionSources
												);*/
	
	return 0;
 
}
