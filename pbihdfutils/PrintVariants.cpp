#include "data/hdf/HDFCmpFile.h"
#include "datastructures/alignment/CmpFile.h"
#include "CommandLineParser.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFBasWriter.h"
#include "utils/RegionUtils.h"
#include "files/ReaderAgglomerate.h"
#include "utils/FileOfFileNames.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "Types.h"
#include "datastructures/matrix/Matrix.h"
#include <string>
#include <ext/hash_map>
#include <map>

#include "InsertedString.h"
#include "Profile.h"
#include "MSAUtils.h"

using namespace std;

typedef vector<unsigned char> ByteAlignment;



int PrintSNVS(FASTASequence &genome, DNASequence &sample, vector<int> &coverage) {
  int i;
  for (i = 0; i < genome.size(); i++) {
    //#ref000001	.	SNV	9454	9454	.	.	.	coverage=40;confidence=10.0;genotype=A;reference=C
    if (genome.seq[i] != sample.seq[i] and sample.seq[i] != ' ') {
      cout << genome.GetName() << "\t.\tSNV\t" << i + 1 << "\t" << i + 1
           << "\t.\t.\t.\tcoverage=" << coverage[i] << ";confidence=10.0;"
           << "variantseq="<<sample.seq[i] <<";"
           << "reference="<<genome.seq[i] << endl;
    }
  }
}

void PrintDeletions(FASTASequence &reference, DNASequence &sample, vector<int> &coverage) {
  //ref000001	.	deletion	50141	50141	.	.	.	coverage=28;confidence=1.0;length=1;reference=C
  int ndel = 0;
  int i;
  for (i = 0; i < reference.length; i++) {
    if (sample.seq[i] == ' ') {
      cout << reference.GetName() << "\t.\tdeletion\t" << i +1 << "\t" << i+1 << "\t.\t.\t.\t"
           << "coverage="<<coverage[i+1]<<";"
           << "confidence=10.0;length=1;reference=" << reference.seq[i] << endl;
    }
  }
}



int CountDeletions(DNASequence &seq) {
  int ndel = 0;
  int i;
  for (i = 0; i < seq.length; i++) {
    if (seq.seq[i] == ' ') {
      ndel++;
    }
  }
  return ndel;
}
  
void StoreCalledBasesInGenome(Profile &profile, DNASequence &genome) {
  int i;
  for (i = 0; i < profile.profileMatrix.GetNCols(); i++) {
    int c;
    int maxC = 0;
    int maxCCount = profile.profileMatrix[0][i];
    for (c = 1; c < 5; c++) {
      if (profile.profileMatrix[c][i] > maxCCount) {
        maxCCount = profile.profileMatrix[c][i];
        maxC = c;
      }
    }
    if (maxC == 4) {
      //      cout << " max del at " << i << " with " << maxCCount << " " << profile.profileMatrix[0][i] << " " << profile.profileMatrix[1][i] << " " << profile.profileMatrix[2][i] << " " << profile.profileMatrix[3][i] << " " << endl;
      genome.seq[i] = ' ';
    }
    else {
      assert(maxC < 4);
      genome.seq[i] = TwoBitToAscii[maxC];
    }
  }
}

int PrintInsertions(Profile &profile, FASTASequence &sampleGenome, vector<int> &coverage, int window=5) {
  //
  // Look to see if any insertion columns appear about average in
  // relation to the average coverage in the region.
  //
  int i;
  int numInsertions=0;
  for (i = 0; i < profile.profileMatrix.GetNCols(); i++) {
    //
    // First find the character inserted the most.
    //
    int maxInsChar = 5;
    int maxInsCount = profile.profileMatrix[5][i];
    int c;
    for (c = 6; c < 9; c++) {
      if (profile.profileMatrix[c][i] > maxInsCount) {
        maxInsCount = profile.profileMatrix[c][i];
        maxInsChar  = c;
      }
    }
    //
    // Now count the coverage in the region.
    //
    int start, end;
    if (i < window) {
      start = 0;
    }
    else {
      start = i - window;
    }
    if (i + window > profile.profileMatrix.GetNCols()) {
      end = profile.profileMatrix.GetNCols();
    }
    else {
      end = i + window;
    }
    int totalCov = 0;
    int w;
    for (w = start; w < end; w++) {
      if (sampleGenome[w] != ' ') {
        totalCov += profile.profileMatrix[TwoBit[sampleGenome[w]]][w];
      }
      else {
        totalCov += profile.profileMatrix[4][w];
      }
    }
    if (end > start) {
      float averageCoverage = totalCov / ((float)(end - start));
      if (maxInsCount * 2 > averageCoverage) {
        //ref000001	.	insertion	50223	50223	.	.	.	coverage=34;confidence=10.0;length=1;variantseq=A;reference=G
        cout << sampleGenome.GetName() << "\t.\tinsertion\t" << i + 1 << "\t" << i + 1 << "\t.\t.\t.\t"
             << "coverage=" << coverage[i] << ";confidence=10.0;length=1;" 
             << "variantseq=" << TwoBitToAscii[maxInsChar-5] << ";"
             << "reference=" << sampleGenome.seq[i] << endl;
        ++numInsertions;
      }
    }                                          
  }
  return numInsertions;
}

void ProfileToSequence(Profile &profile, DNASequence &seq) {




}


typedef vector<InsertedString> InsertedStringList;
void StoreInsertionQVs(InsertedStringList &insertions, vector<char> &valueStr) {
  int i;
  for (i = 0; i < insertions.size(); i++) {
    int substrLen = insertions[i].insSeq.size();
    int substrPos = insertions[i].alnPos;
    // Modify the string to be from the valuestr.
    insertions[i].insQVSeq = "";
    insertions[i].insQVSeq.insert(0, &valueStr[substrPos], substrLen);
  }
}



int ComputeOffsetToRefPosition(ByteAlignment &aln, DNALength alnStartPos, DNALength refPos, DNALength &alnIndex, DNALength &refOffset) {
  alnIndex = 0;
  refOffset = 0;
  while(alnIndex < aln.size() and refOffset + alnStartPos != refPos) {
    // move past the current position.
    if (RefChar[aln[alnIndex]] != ' ') {
      refOffset++;
    }
    alnIndex++;
  }
  if (alnIndex < aln.size()) {
    return 1;
  }
  else {
    return 0;
  }
}


int main(int argc, char* argv[]) {
  CommandLineParser clp;
  string cmpH5FileName, genomeFileName;
  int  scale = 5;
  int  marginSize = 5;
  stringstream verboseHelpStream, helpStream;
  vector<int> coverage;
  verboseHelpStream << "printmsa is a utility to either view the msa centered at a position" << endl
                    << "in a base.h5 file, or print all reads overlapping a position." << endl;
  clp.SetVerboseHelp(verboseHelpStream.str());
  clp.RegisterStringOption("cmpH5FileName", &cmpH5FileName, "The input hdf file");
  clp.RegisterStringOption("genome", &genomeFileName, "Name of the reference reads were aligned to (1 contig for now)");
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterIntOption("scale",           &scale, "Scale quality values by 'scale'", CommandLineParser::PositiveInteger);

  clp.ParseCommandLine(argc, argv);

  FASTASequence genome;
  ReaderAgglomerate reader;
  reader.Initialize(genomeFileName);

  
  CmpFile cmpFile;
	HDFCmpFile<CmpAlignment> cmpReader;


	if (cmpReader.Initialize(cmpH5FileName, H5F_ACC_RDONLY) == 0) {
		cout << "ERROR, could not open the cmp file." << endl;
		exit(1);
	}
  //	cout << "Reading cmp file." << endl;
	cmpReader.Read(cmpFile);

  int p;
  while (reader.GetNext(genome)) {
    genome.ToUpper();
    FASTASequence sampleGenome;
    sampleGenome.Copy(genome);
    Profile profile(genome.length);
    coverage.resize(genome.length);
    fill(coverage.begin(), coverage.end(), 0);
    int alignmentIndex;
    
    for (alignmentIndex = 0; alignmentIndex < cmpFile.alnInfo.alignments.size(); alignmentIndex++) {
      int readStart, readEnd, refStart, refEnd;

      readStart = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart();
      readEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd();
      refStart  = cmpFile.alnInfo.alignments[alignmentIndex].GetRefStart();
      refEnd    = cmpFile.alnInfo.alignments[alignmentIndex].GetRefEnd();
      int holeNumber = cmpFile.alnInfo.alignments[alignmentIndex].GetHoleNumber();
      int refGroupId = cmpFile.alnInfo.alignments[alignmentIndex].GetRefGroupId();
      int refInfoIndex;
      int refGroupIndex = cmpReader.refGroupIdToArrayIndex[refGroupId];
      int refInfoId = cmpFile.refGroup.refInfoId[refGroupIndex];
      if (cmpFile.refInfo.RefIdToIndex(refInfoId, refInfoIndex) == false) {
        cout << "ERROR, could not find refGroupId " << refGroupId << " in ref info structure." << endl;
        exit(1);
      }
      
      string refFullName = cmpFile.refInfo.fullName[refInfoIndex];
      //      cout <<" The alignment is aligned to " << refFullName << " vs " << genome.title << endl;
      
      int refStrand      = cmpFile.alnInfo.alignments[alignmentIndex].GetRCRefStrand();

      //
      // Read the alignment string.  All alignments 
      //
      int offsetBegin = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
      int offsetEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd();
		
      int alignedSequenceLength = offsetEnd - offsetBegin;
      string   alignedSequence;
      ByteAlignment byteAlignment, byteAlignmentRC;
      vector<char> insQVChars, delQVChars, subQVChars;
      
      if (refFullName != genome.title) {
        //  cout << "ref full name is not " << genome.title << endl;
        continue;
      }
      if (alignedSequenceLength >= 0) {
        alignedSequence.resize(alignedSequenceLength);
        byteAlignment.resize(alignedSequenceLength);
      }
      byteAlignment = cmpFile.alnInfo.alignments[alignmentIndex].alignmentArray;
    
      if (refStrand == 1) {
        byteAlignmentRC.resize(byteAlignment.size());
        MakeReverseComplementByteAlignment(&byteAlignment[0], byteAlignment.size(), &byteAlignmentRC[0]);
        byteAlignment = byteAlignmentRC;
      }

      vector<int> refPositions, queryPositions;
      ComputeQueryPositions(byteAlignment, queryPositions);
      ComputeRefPositions(byteAlignment,   refPositions);
    
      string alnStr;
      InsertedStringList insertions;
      StoreInsertedStrings(byteAlignment, refPositions, queryPositions,
                           insertions, 0, byteAlignment.size());
            
      //
      // Project the alignment onto the reference
      int p = 0, bp = 0;
      while (bp < byteAlignment.size()) {
        if (RefChar[byteAlignment[bp]] != ' ') {
          alnStr.push_back(QueryChar[byteAlignment[bp]]);
          p++;
        }
        bp++;
      }

      for (p = 0; p < refPositions.size(); p++) {
        if (refPositions[p] != -1) {
          coverage[refStart + refPositions[p]]++;
        }
      }

      profile.StoreProfile(alnStr, 0, 0, insertions, refStart);
      /*    if (alignmentIndex % 10000 == 0) {
            cout << alignmentIndex << endl;
            }*/
    }
    StoreCalledBasesInGenome(profile, sampleGenome);
    PrintDeletions(genome, sampleGenome, coverage);
    PrintSNVS(genome, sampleGenome, coverage);
    PrintInsertions(profile, sampleGenome, coverage);
  }
  //  cout << "The consensus has " << ndel << " deletions. " << endl;
//  profile.Print(6);
  return 0;
}
