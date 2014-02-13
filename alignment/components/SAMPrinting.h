#ifndef ALIGNMENT_COMPONENTS_SAM_PRINTING_H_
#define ALIGNMENT_COMPONENTS_SAM_PRINTING_H_

#include <vector>

#include "utils/MD5Utils.h"
#include "utils/StringUtils.h"
#include "datastructures/metagenome/SequenceIndexDatabase.h"

void MakeSAMHDString(string majorVersion, string &hdString) {
  stringstream out;
  out << "@HD\t" << "VN:" << majorVersion;
  hdString = out.str();
}

void MakeSAMPGString(string commandName, string version, 
                     string &commandLineString, string &pgString) {
  stringstream out;
  out << "@PG\t" << "ID:" << commandName << "\tVN:"<< version <<"\tCL:"<<commandLineString;
  pgString = out.str();
}

void ParseChipIdFromMovieName(string &movieName, string &chipId) {
  vector<string> movieNameFields;
  Tokenize(movieName, "_", movieNameFields);
  if (movieNameFields.size() == 1) {
    chipId = movieNameFields[0];
  }
  else if (movieNameFields.size() > 4) {
    chipId = movieNameFields[3];
  }
  else {
    chipId = "NO_CHIP_ID";
  }
}


void PrepareSAMHeader(string programName, string majorVersion, string version,
                      ReaderAgglomerate *reader,
                      vector<string> &readsFileNames, 
                      SequenceIndexDatabase<FASTASequence> &seqdb, 
                      string commandLineString, 
                      ostream &outFile) {
  string hdString, sqString, rgString, pgString;
  MakeSAMHDString(majorVersion, hdString);
  outFile << hdString << endl;
  seqdb.MakeSAMSQString(sqString);
  outFile << sqString; // this already outputs endl
  
  //
  //  Look through the fofn to list all of the read group names by
  //  querying the bas.h5 files.  Files that are not bas.h5 use the
  //  file name.
  //
  int readsFileIndex;
  for (readsFileIndex = 0; readsFileIndex < readsFileNames.size(); readsFileIndex++ ) {    
    reader->SetReadFileName(readsFileNames[readsFileIndex]);
    reader->Initialize();
    string movieName, movieNameMD5;
    reader->GetMovieName(movieName);
    MakeMD5(movieName, movieNameMD5, 10);
    string chipId;
    ParseChipIdFromMovieName(movieName, chipId);
    // 
    // If this isn't a real movie name, the chipid is a placeholder,
    // and is not meaningful.
    //
    outFile << "@RG\t" << "ID:"<<movieNameMD5<<"\t" << "PU:"<< movieName << "\tSM:"<<chipId << endl;
    reader->Close();
  }

  //
  // Record what the command line was to do the alignments.
  //
  MakeSAMPGString(programName, version, commandLineString, pgString);
  outFile << pgString << endl;
}


#endif
