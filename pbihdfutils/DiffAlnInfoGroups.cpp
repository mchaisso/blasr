#include "data/hdf/HDFCmpFile.h"
#include "datastructures/alignment/CmpFile.h"
#include "CommandLineParser.h"

void CopyPointers(vector<CmpAlignment> &alignments, vector<CmpAlignment*> &alnPtrs) {
  alnPtrs.resize(alignments.size());
  int i;
  for (i = 0; i < alignments.size(); i++) {
    alnPtrs[i] = &alignments[i];
  }
}

class CompareAlnPtrs {
public:
  int operator()(CmpAlignment *aPtr, CmpAlignment *bPtr) {
    return *aPtr < *bPtr;
  }
};
  
int main(int argc, char* argv[]) {

  CommandLineParser clp;
  string cmpFileAName, cmpFileBName;
  int level = 0;
  bool noSort = false;
  clp.RegisterStringOption("cmpA", &cmpFileAName, "First cmp file to compare.");
  clp.RegisterStringOption("cmpB", &cmpFileBName, "Second cmp file to compare.");
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterIntOption("level", &level, "How deep to align, 1=count aligned reads"
                        " 2 = count aligned bases, 3 = compare boundaries, 4 = "
                        " compare alignment strings.", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("nosort", &noSort, "Do not sort the cmpH5 files.");

  clp.ParseCommandLine(argc, argv);

  HDFCmpFile<CmpAlignment> cmpReader;
  CmpFile cmpFileA, cmpFileB;

  //
  // Read in the two cmp files. Not for large datasets.
  //
	if (cmpReader.Initialize(cmpFileAName, H5F_ACC_RDONLY) == 0) {
    cout << "Could not open " << cmpFileAName << endl;
  }

  cmpReader.Read(cmpFileA);

  cmpReader.Close();

	if (cmpReader.Initialize(cmpFileBName, H5F_ACC_RDONLY) == 0) {
    cout << "Could not open " << cmpFileBName << endl;
  }
  cmpReader.Read(cmpFileB);

  if (level >= 0) {
    if (cmpFileA.alnInfo.alignments.size() == cmpFileB.alnInfo.alignments.size()) {
      cout << "LEVEL 0 same number of alignments" << endl;
    }
    else {
      cout << "LEVEL 0 " << cmpFileA.alnInfo.alignments.size() << " " << cmpFileB.alnInfo.alignments.size() << endl;
      //      return 0;
    }
    if (level == 0) { return 0; }
  }

  int i;

  if (level >= 1) {
    int nABases = 0, nBBases = 0;
    for (i = 0; i < cmpFileA.alnInfo.alignments[i].alignmentIndex.size(); i++) {
      nABases += cmpFileA.alnInfo.alignments[i].alignmentIndex[12] - cmpFileA.alnInfo.alignments[i].alignmentIndex[11];
    }
    for (i = 0; i < cmpFileB.alnInfo.alignments[i].alignmentIndex.size(); i++) {
      nBBases += cmpFileB.alnInfo.alignments[i].alignmentIndex[12] - cmpFileB.alnInfo.alignments[i].alignmentIndex[11];
    }
    if (nABases == nBBases) {
      cout << "LEVEL 1 same number of bases." << endl;
    }
    else {
      cout << "LEVEL 1 different number of bases " << nABases << " " << nBBases << "." << endl;
    }
    if (level == 1) { return 0; }
  }
  vector<CmpAlignment*> cmpAlignmentsA, cmpAlignmentsB;

  CopyPointers(cmpFileA.alnInfo.alignments, cmpAlignmentsA);
  CopyPointers(cmpFileB.alnInfo.alignments, cmpAlignmentsB);
  if (noSort == false) {
    cout << "sorting a" << endl;
    sort(cmpAlignmentsA.begin(), cmpAlignmentsA.end(), CompareAlnPtrs());
    cout << "sorting b" << endl;
    sort(cmpAlignmentsB.begin(), cmpAlignmentsB.end(), CompareAlnPtrs());
  }
  int a;
  if (level >= 2) {
    for (a = 0; a < cmpAlignmentsA.size() and a < cmpAlignmentsB.size(); a++) {
      if ((*cmpAlignmentsA[a]).alignmentIndex[7] != (*cmpAlignmentsB[a]).alignmentIndex[7] or
          (*cmpAlignmentsA[a]).alignmentIndex[4] != (*cmpAlignmentsB[a]).alignmentIndex[4] or
          (*cmpAlignmentsA[a]).alignmentIndex[5] != (*cmpAlignmentsB[a]).alignmentIndex[5] or
          (*cmpAlignmentsA[a]).alignmentIndex[11] != (*cmpAlignmentsB[a]).alignmentIndex[11] or
          (*cmpAlignmentsA[a]).alignmentIndex[12] != (*cmpAlignmentsB[a]).alignmentIndex[12]) {
        cout << "LEVEL 2 is different, nuff said for now." << endl;
        return 0;
      }
    }
    cout << "LEVEL 2 same coordinates for each sorted alignment." << endl;
  }

  if (level >= 3) {
    for (a = 0; a < cmpAlignmentsA.size() and a < cmpAlignmentsB.size(); a++) {
      if (memcmp(&(*cmpAlignmentsA[a]).alignmentArray[0], &(*cmpAlignmentsB[a]).alignmentArray[0],  (*cmpAlignmentsB[a]).alignmentArray.size()) != 0) {
        cout << "LEVEL 3 is different, nuff said for now." << endl;
        return 0;
      }
    }
    cout << "LEVEL 3 same alignment arrays." << endl;
  }
 
  return 1;
}
