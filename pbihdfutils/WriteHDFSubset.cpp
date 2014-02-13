#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "data/hdf/HDFBasWriter.h"
#include "data/hdf/PlatformId.h"
#include "utils/StringUtils.h"
#include "utils/RegionUtils.h"
#include "utils/FileOfFileNames.h"
#include "SMRTSequence.h"
#include <string>
#include <iostream>
#include <vector>

void PrintUsage() {
		cout << "usage: writeHDFSubset in out [idx1 idx2 idx3] [-pat p] [-fromto from to] [-table] [-v] ..." << endl;
		cout << " idx       is the index of a read (starting at 0), " << endl
				 << "             no indices needed when in is a file index (.findex) file." << endl
				 << " pat p     is a pattern to extract all reads with p in their title." << endl
         << " fromto    is a range where from < to." << endl
				 << " If 'in' is a file index,  " << endl;
}


int main(int argc, char* argv[]) {
	
	string inFileName, outFileName;

	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	inFileName  = argv[1];
	outFileName = argv[2];

	vector<int> readIndices;
	int argi = 3;
	vector<string> patterns;
	vector<int> holeNumbers;
  string regionTableFileName = "";
  int from = 0, to = 0;

	map<string, vector<int> > inputFiles;
	bool useFileIndex = false;
	vector<string> inFiles;
	if (inFileName.find("fofn") != inFileName.npos) {
		//
		// Input file is a fofn, read in file names.
		//
		FileOfFileNames::StoreFileOrFileList(inFileName, inFiles);
		int i;
		for (i = 0; i < inFiles.size(); i++) {
			inputFiles[inFiles[i]] = vector<int>(0);
		}
	}
	else if (inFileName.find("findex") != inFileName.npos) {
		useFileIndex = true;
		ifstream fileIndex;
		CrucialOpen(inFileName, fileIndex, std::ios::in);
		while (fileIndex) {
			string fileName;
			int index;
			if (!(fileIndex >> fileName >> index)){ 
				break;
			}
			if (fileName == "") {
				break;
			}
			inputFiles[fileName].push_back(index);
		}
		map<string, vector<int> >::iterator mapIt, mapEnd;
		for (mapIt = inputFiles.begin(); mapIt != inputFiles.end(); ++mapIt) {
			std::sort((*mapIt).second.begin(), (*mapIt).second.end());
			inFiles.push_back((*mapIt).first);
		}
	}
	else {
		inputFiles[inFileName] = vector<int>(0);
	}

	bool verbose = false;

	while (argi < argc) {
		if (strlen(argv[argi]) > 0 and argv[argi][0] == '-'){ 
			if (strcmp(argv[argi], "-holeNumber") == 0) {
				holeNumbers.push_back(atoi(argv[++argi]));
			}
			else if (strcmp(argv[argi], "-regionTable") == 0) {
				regionTableFileName = argv[++argi];
			}
      else if (strcmp(argv[argi], "-fromto") == 0) {
        from = atoi(argv[++argi]);
        to   = atoi(argv[++argi]);
        if (from >= to) {
          cout <<"ERROR. From must be less than to." << endl;
          exit(1);
        }
      }
			else if (strcmp(argv[argi], "-v") == 0) {
				verbose = true;
			}
      else {
        cout <<"ERROR. Bad option " << argv[argi] << endl;
        PrintUsage();
        exit(1);
      }
		}
		else {
			readIndices.push_back(atoi(argv[argi]));
		}
		++argi;
	}

	//
	// For many possible reads from one input file.
	//
  int index;
  for (index = from; index < to; index++) {
    readIndices.push_back(index);
  }
	std::sort(readIndices.begin(), readIndices.end());
	map<string, vector<int> >::iterator mapIt;
	if (readIndices.size() > 0) {
		for (mapIt = inputFiles.begin(); mapIt != inputFiles.end(); ++mapIt) {
			(*mapIt).second = readIndices;
		}
	}

	T_HDFBasReader<SMRTSequence> reader;
  HDFRegionTableReader regionReader;

	HDFBasWriter writer;
  HDFRegionTableWriter regionWriter;
	reader.InitializeDefaultIncludedFields();
	writer.InitializeDefaultIncludedFields();
	writer.IncludeField("HoleNumber");
	writer.IncludeField("HoleXY");

	bool outputIsInitialized = false;

		
	int inFileIndex;

	RegionTable regionTable;

	for (mapIt = inputFiles.begin(); mapIt != inputFiles.end(); ++mapIt) {
		inFileName = (*mapIt).first;

		reader.Initialize(inFileName);
		regionReader.Initialize(inFileName);
		regionReader.ReadTable(regionTable);

		string changeListID;
		reader.GetChangeListID(changeListID);

		if (outputIsInitialized == false) {
			if (reader.scanDataReader.GetPlatformId() == AstroPlatform) {
				writer.Initialize(outFileName, reader.GetMovieName(), reader.GetRunCode());
			}
			else {
				writer.Initialize(outFileName, reader.GetMovieName(), changeListID);
			}
			regionWriter.Initialize(writer.pulseDataGroup);
			outputIsInitialized = true;
		}


		int ri;
		int curReadIndex = 0;
		SMRTSequence seq;
		bool printSeq = false;
		ri = 0;
		reader.PrepareForRandomAccess();
		readIndices = (*mapIt).second;
		for (ri = 0; ri < readIndices.size(); ri++) {
			reader.GetReadAt(readIndices[ri], seq);
			if (verbose) {
				cout << seq.title << endl;
			}
			writer.Write(seq);

			//
			// Write out region information for the read.
			//
			int low, high;
			FindRegionIndices(readIndices[ri], &regionTable, low, high);
			int regionIndex;
			for (regionIndex = low; regionIndex < high; regionIndex++) {
				regionWriter.Write(regionTable.table[regionIndex]);
			}
		}
  }

	regionWriter.Finalize(regionTable.columnNames,
												regionTable.regionTypes, 
												regionTable.regionDescriptions, 
												regionTable.regionSources
												);

	writer.Flush();
}
