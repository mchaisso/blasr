#include "data/hdf/HDFCmpFile.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "FASTASequence.h"
#include <string>


using namespace std;

int main() {
  string cmpFileName = "test.cmp.h5";
  
  HDFCmpFile<AlignmentCandidate<FASTASequence, FASTASequence> > cmpFile;

  cmpFile.Create(cmpFileName);
  cmpFile.alnGroupGroup.AddPath("/ref000001/m111217_220951_richard_c100227132550000001523001504251270_s1_p0", 1);
  cmpFile.alnGroupGroup.AddPath("/ref000001/some_other_movie", 2);
}
  
