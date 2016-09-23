#include "datastructures/alignment/AlignmentCandidate.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "data/hdf/HDFCmpFile.h"
#include "FASTASequence.h"
#include "FASTAReader.h"
#include "CommandLineParser.h"
#include "datastructures/alignmentset/AlignmentSetToCmpH5Adapter.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "utils/ChangeListID.h"
#include "utils/TimeUtils.h"
#include <iostream>

char VERSION[] = "v1.0.0";
char PERFORCE_VERSION_STRING[] = "$Change: 126414 $";

int main(int argc, char* argv[]) {
  string program = "samtoh5";
  string versionString = VERSION;
  AppendPerforceChangelist(PERFORCE_VERSION_STRING, versionString);
  string samFileName, cmpFileName, refFileName;
  bool parseSmrtTitle = false;
  bool useShortRefName = false;
  CommandLineParser clp;
  string readType = "standard";
  int verbosity = 0;

  clp.SetProgramName(program);
  clp.SetProgramSummary("Converts in.sam file to out.cmp.h5 file.");
  clp.SetVersion(versionString);

  clp.RegisterStringOption("in.sam", &samFileName, 
                           "Input SAM file.", true);
  clp.RegisterStringOption("reference.fasta", &refFileName, 
                           "Reference used to generate reads.", true);
  clp.RegisterStringOption("out.cmp.h5", &cmpFileName, 
                           "Output cmp.h5 file.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterFlagOption("smrtTitle", &parseSmrtTitle, 
                         "Use this option when converting alignments "
                         "generated from reads produced by the "
                         "pls2fasta from bas.h5 files by parsing read "
                         "coordinates from the SMRT read title.  The title " 
                         "is in the format /name/hole/coordinates, where "
                         "coordinates are in the format \\d+_\\d+, and "
                         "represent the interval of the read that was "
                         "aligned.");
  clp.RegisterStringOption("readType", &readType, 
                         "Set the read type: 'standard', 'strobe', 'CCS', "
                         "or 'cDNA'");
  clp.RegisterIntOption("verbosity", &verbosity, 
                         "Set desired verbosity.", 
                         CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("useShortRefName", &useShortRefName, 
                         "Use abbreviated reference names obtained "
                         "from file.sam instead of using full names "
                         "from reference.fasta.");
  string description = ("Because SAM has optional tags that have different "
    "meanings in different programs, careful usage is required in order to "
    "have proper output. The \"xs\" tag in bwa-sw is used to show the "
    "suboptimal score, but in PacBio SAM (blasr) it is defined as the start "
    "in the query sequence of the alignment.\nWhen \"-smrtTitle\" is "
    "specified, the xs tag is ignored, but when it is not specified, the "
    "coordinates given by the xs and xe tags are used to define the interval "
    "of a read that is aligned. The CIGAR string is relative to this interval.");
  clp.SetExamples(description);

  clp.ParseCommandLine(argc, argv);

  if (readType != "standard" and readType != "strobe" and 
      readType != "cDNA" and readType != "CCS") {
    cout << "ERROR. Read type '" << readType 
         << "' must be one of either 'standard', 'strobe', 'cDNA' or 'CCS'." 
         << endl;
    exit(1);
  }
    
  cerr << "[INFO] " << GetTimestamp() << " [" << program << "] started." << endl;

  SAMReader<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> samReader;
  FASTAReader fastaReader;
  HDFCmpFile<AlignmentCandidate<FASTASequence, FASTASequence> > cmpFile;

  //
  // Initialize input/output files.
  //
  samReader.Initialize(samFileName);
  fastaReader.Initialize(refFileName);
  cmpFile.Create(cmpFileName);

  //
  // Configure the file log.
  //
  string command;
  CommandLineParser::CommandLineToString(argc, argv, command);
  string log = "Convert sam to cmp.h5";
  cmpFile.fileLogGroup.AddEntry(command, log, program, GetTimestamp(), versionString);

  //
  // Set the readType
  //
  cmpFile.SetReadType(readType);

  //
  // Read necessary input.
  //

  vector<FASTASequence> references;
	fastaReader.storeName = true;
  fastaReader.ReadAllSequences(references);
  
  //
  // This should probably be handled by the alignmentSetAdapter, but
  // time constraints...
  //
  AlignmentSet<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> alignmentSet;
  samReader.ReadHeader(alignmentSet);
 
  //
  // The order of references in vector<FASTASequence> references and
  // AlignmentSet<, , >alignmentSet.references can be different.
  // Rearrange alignmentSet.references such that it is ordered in
  // exactly the same way as vector<FASTASequence> references.
  //
  alignmentSet.RearrangeReferences(references);

  //
  // Always recompute the MD5 values even if they exist in the input
  // sam file. Because MD5 is defined differently in sam and cmp.h5 files.
  // The SAM convention uppercases and normalizes before computing the MD5. 
  // For cmp.h5, we compute the MD5 on the sequence 'as is'.
  // 
  for(int i = 0; i < alignmentSet.references.size(); i++) {
      MakeMD5((const char*)&references[i].seq[0], 
              (unsigned int)references[i].length, alignmentSet.references[i].md5);
  }
 
  //
  // Map short names for references obtained from file.sam to full names obtained from reference.fasta
  //
  map<string, string> shortRefNameToFull;
  map<string, string>::iterator it;
  if (references.size() != alignmentSet.references.size()) {
		cout<< "WARNING. The references in the sam file header do not match those in the supplied reference." << endl;
		
	}
  if (!useShortRefName) {
		map<string, int> referenceNameToIndex;
		for (int i = 0; i < references.size(); i++) {
			referenceNameToIndex[references[i].GetName()] = i;
		}
		for (int i = 0; i < alignmentSet.references.size(); i++) {
			string shortRefName = alignmentSet.references[i].GetSequenceName();
			string fullRefName(references[referenceNameToIndex[shortRefName]].title); 
			if (shortRefNameToFull.find(shortRefName) != shortRefNameToFull.end()) {
				cout << "ERROR, Found more than one reference " << shortRefName << "in sam header" << endl;
				exit(1);
			} 
			shortRefNameToFull[shortRefName] = fullRefName;
			alignmentSet.references[i].sequenceName = fullRefName;
		}
  }

  //
  // Start setting up the cmp.h5 file.
  //
  AlignmentSetToCmpH5Adapter<HDFCmpFile<AlignmentCandidate<FASTASequence, FASTASequence> > > alignmentSetAdapter;
  alignmentSetAdapter.Initialize(references);
	
	alignmentSetAdapter.StoreAllMovieInfo(alignmentSet.readGroups, cmpFile);
  
  //
  // Store the alignments.
  //
  SAMAlignment samAlignment;
  int alignIndex = 0;
	
  while (samReader.GetNextAlignment(samAlignment)) {
    if (samAlignment.rName == "*") {
      continue;
    }
    if (!useShortRefName) {
        //convert shortRefName to fullRefName
        it = shortRefNameToFull.find(samAlignment.rName);
        if (it == shortRefNameToFull.end()) {
            cout << "ERROR, Could not find " << samAlignment.rName << " in the reference repository." << endl;
            exit(1);
        }
        samAlignment.rName = (*it).second;
    }
    vector<AlignmentCandidate<> > convertedAlignments;
    if (verbosity > 0) {
      cout << "Storing alignment for " << samAlignment.qName << endl;
    }
		
		if (alignmentSetAdapter.ReferenceIsStored(samAlignment.rName) == false) {
			int refIndex = alignmentSet.refNameToIndex[samAlignment.rName];
			alignmentSetAdapter.StoreReferenceInfo(alignmentSet.references[refIndex], cmpFile);			
		}
    SAMAlignmentsToCandidates(samAlignment, 
                              references, alignmentSetAdapter.refNameToAllRefIndex,
                              convertedAlignments, parseSmrtTitle, false);
		

    alignmentSetAdapter.StoreAlignmentCandidateList(convertedAlignments, cmpFile, alignIndex);
    int a;
    for (a = 0; a < convertedAlignments.size(); a++) {
      convertedAlignments[a].FreeSubsequences();
    }
    ++alignIndex;
    /*    if (alignIndex == 100) {
      return 0;
      }*/
  }

  cerr << "[INFO] " << GetTimestamp() << " [" << program << "] ended." << endl;
  return 0;
}
