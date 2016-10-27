#include "datastructures/alignment/AlignmentCandidate.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "data/hdf/HDFBasWriter.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "CommandLineParser.h"
#include "datastructures/alignmentset/AlignmentSetToCmpH5Adapter.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "datastructures/alignmentset/SAMAlignment.h"
#include "datastructures/alignmentset/SAMToSMRTSequence.h"
#include "datastructures/alignmentset/SAMQVConversion.h"
#include "SMRTSequence.h"
#include "datastructures/reads/RegionTable.h"
#include "datastructures/reads/BaseFile.h"
#include "utils/ChangeListID.h"
#include "utils/TimeUtils.h"
#include <iostream>



int main(int argc, char* argv[]) {
  string program = "samtobas";


  string samFileName, basFileName;
  CommandLineParser clp;

  int verbosity = 0;
	bool defaultToP6 = false;
  clp.SetProgramName(program);
  clp.SetProgramSummary("Converts in.sam file to out.cmp.h5 file.");


  clp.RegisterStringOption("in.sam", &samFileName, 
                           "Input SAM file.  Use /dev/stdin for piping.", true);
  clp.RegisterStringOption("out.bas.h5", &basFileName, 
                           "Output bas.h5 file.", true);
	
  clp.RegisterPreviousFlagsAsHidden();
	clp.RegisterFlagOption("defaultToP6", &defaultToP6, "Default to P6 when changelist,chemistry,kit triple is not present", false);
  string description = ("Transform data in a sam file into a bas.h5 file. "
												" This will create fake hole numbers and movie names for the reads, but"
												" they should otherwise be usable by any downstream pipeline as if they "
												" were native bas.h5 files from a machine, minus IPD values.");

  clp.SetExamples(description);

  clp.ParseCommandLine(argc, argv);

	cerr << "[INFO] " << GetTimestamp() << " [" << program << "] started." << endl;

  SAMReader<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> samReader;
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


  //
  //
  AlignmentSet<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> alignmentSet;
  samReader.ReadHeader(alignmentSet);
	string changelistId, bindingKit, sequencingKit;
	alignmentSet.GetRepresentativeChangelistId(changelistId);
	alignmentSet.GetRepresentativeBindingKit(bindingKit);
	alignmentSet.GetRepresentativeSequencingKit(sequencingKit);
	
	if (defaultToP6 and (bindingKit == "" or changelistId == "" or sequencingKit == "")) {
		cerr << "WARNING, providing default run data for P6 chemistry." << endl;
		//
		// Provide default
		bindingKit = "100356300";
		sequencingKit = "100356200";
		changelistId = "2.3.0.0.140640";
	}	

	writer.Initialize(basFileName, "pileup_reads", changelistId);
	writer.AddBindingKit(bindingKit);
	writer.AddSequencingKit(sequencingKit);
	regionWriter.Initialize(writer.pulseDataGroup);
  
  //
  // Convert alignments to reads.  Separate subread information is discarded.
  //
  SAMAlignment samAlignment;
	int alignmentIndex = 0;

  while (samReader.GetNextAlignment(samAlignment)) {
    if (samAlignment.rName == "*") {
      continue;
    }
		samAlignment.zmw = alignmentIndex;		

		SMRTSequence read;
		ConvertSAMToSMRTSequence(samAlignment, read);		

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

