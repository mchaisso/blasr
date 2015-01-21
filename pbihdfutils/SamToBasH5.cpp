#include "datastructures/alignment/AlignmentCandidate.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "data/hdf/HDFBasWriter.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "CommandLineParser.h"
#include "datastructures/alignmentset/AlignmentSetToCmpH5Adapter.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "datastructures/alignmentset/SAMQVConversion.h"
#include "datastructures/reads/RegionTable.h"
#include "datastructures/reads/BaseFile.h"
#include "utils/ChangeListID.h"
#include "utils/TimeUtils.h"
#include <iostream>


int CopyFieldValues(unsigned char* dest, string src) {
	if (src.size() == 0) {
		return 0;
	}
	else {
		memcpy(dest, src.c_str(), src.size());
		return src.size();
	}
}



int AssignData( unsigned char* &field, string src) {
	if (field == NULL) {
		if (src.size() > 0) {
			cout << "ERROR. Copying into unallocated field." << endl;
			assert(0);
		}
		else {
			return 0;
		}
	}
	int res = CopyFieldValues(field, src);
	return res;
}

int AssignQualityData( unsigned char* &field, string src) {
	int res = AssignData(field, src);

	if (res == 0) {
		field = NULL;
	}		
	else {
		QualityStringToStored(field, src.size());
	}
	return res;
}

int main(int argc, char* argv[]) {
  string program = "samtobas";


  string samFileName, basFileName;
  CommandLineParser clp;

  int verbosity = 0;

  clp.SetProgramName(program);
  clp.SetProgramSummary("Converts in.sam file to out.cmp.h5 file.");


  clp.RegisterStringOption("in.sam", &samFileName, 
                           "Input SAM file.  Use /dev/stdin for piping.", true);
  clp.RegisterStringOption("out.bas.h5", &basFileName, 
                           "Output bas.h5 file.", true);
  clp.RegisterPreviousFlagsAsHidden();
  string description = ("Transform data in a sam file into a bas.h5 file. "
												" This will create fake hole numbers and movie names for the reads, but"
												" they should otherwise be usable by any downstream pipeline as if they "
												" were native bas.h5 files from a machine, minus IPD values.");

  clp.SetExamples(description);

  clp.ParseCommandLine(argc, argv);

	//  cerr << "[INFO] " << GetTimestamp() << " [" << program << "] started." << endl;

  SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> samReader;
	HDFBasReader basReader;
	HDFBasWriter writer;
  HDFRegionTableWriter regionWriter;
  //
  // Initialize input/output files.
  //
  samReader.Initialize(samFileName);


	writer.IncludeField("Basecall");
	writer.IncludeField("DeletionQV");
	writer.IncludeField("DeletionTag");
	writer.IncludeField("InsertionQV");
	writer.IncludeField("SubstitutionTag");
	writer.IncludeField("SubstitutionQV");
	writer.IncludeField("QualityValue");
	writer.IncludeField("HoleNumber");    
	writer.IncludeField("HoleStatus");
	writer.IncludeField("MergeQV");
	writer.IncludeField("HoleNumber"); 
	writer.IncludeField("HoleXY");

	writer.Initialize(basFileName, "pileup_reads", "0.1-beta");
	regionWriter.Initialize(writer.pulseDataGroup);

  //
  // This is not needed for bas.h5 writing, but in order to get the
  // file pointer to the right position, read the header.
  //
  AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet;
  samReader.ReadHeader(alignmentSet);
  
  //
  // Convert alignments to reads.  Separate subread information is discarded.
  //
  SAMAlignment samAlignment;
	int alignmentIndex = 0;
  while (samReader.GetNextAlignment(samAlignment)) {
    if (samAlignment.rName == "*") {
      continue;
    }
		//
		// Populate sequence / quality information.  Any empty qualities
		// are skipped.
		//
		SMRTSequence read;
		read.Allocate(samAlignment.seq.size());
		read.length = samAlignment.seq.size();
		AssignData(read.seq, samAlignment.seq);
		AssignQualityData(read.qual.data, samAlignment.qual);
		AssignQualityData(read.deletionQV.data, samAlignment.qd);
		AssignQualityData(read.insertionQV.data, samAlignment.qi);
		AssignQualityData(read.substitutionQV.data, samAlignment.qs);
		AssignQualityData(read.mergeQV.data, samAlignment.qm);
		AssignData((unsigned char*&) read.substitutionTag, samAlignment.ts);
		AssignData((unsigned char*&) read.deletionTag, samAlignment.td);
		
		read.StoreHoleNumber(alignmentIndex);
		read.StoreHoleStatus(0);
		read.zmwData.numEvents = read.length;


		writer.Write(read);
		RegionAnnotation region;
		region.row[0] = alignmentIndex;
		region.row[1] = 1;
		region.row[2] = 0;
		region.row[3] = read.length;
		region.row[4] = -1;
		// Write the coordinates of the insert.
		regionWriter.Write(region);
		region.row[1] = 2;
		region.row[4] = 1000;
		// Write the high-quality reigon, the same as the insert.
		regionWriter.Write(region);
		
    ++alignmentIndex;
  }
	RegionTable regionTable;
  regionTable.CreateDefaultAttributes();
	regionWriter.Finalize(regionTable.columnNames,
												regionTable.regionTypes, 
												regionTable.regionDescriptions, 
												regionTable.regionSources
												);
	writer.Flush();


  return 0;
}
