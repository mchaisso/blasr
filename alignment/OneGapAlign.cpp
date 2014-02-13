#include "FASTASequence.h"
#include "FASTAReader.h"
#include "algorithms/alignment/OneGapAlignment.h"
#include "CommandLineParser.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "Types.h"
#include "algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "algorithms/alignment/ScoreMatrices.h"

int main(int argc, char* argv[]) {

  string queryFileName, targetFileName;
  if (argc < 3) {
    cout << "usage: oneGapAlign query target " << endl;
    exit(1);
  }

  FASTAReader reader;
  FASTASequence query, target;

  queryFileName = argv[1];
  targetFileName = argv[2];
  reader.Initialize(queryFileName);
  reader.GetNext(query);
  reader.Close();
  reader.Initialize(targetFileName);
  reader.GetNext(target);
  reader.Close();

  FASTASequence leftTarget, rightTarget;
  UInt leftTargetLength = min(target.length, query.length);
  leftTarget.ReferenceSubstring(target, 0, leftTargetLength);

  UInt rightTargetLength = min(target.length - leftTargetLength, query.length);
  rightTarget.ReferenceSubstring(target, target.length - rightTargetLength, rightTargetLength);
  
  DNALength distanceBetweenLeftAndRight = target.length - rightTargetLength - leftTargetLength;

  assert(distanceBetweenLeftAndRight >= 0);

  Alignment alignment;
  alignment.qName = query.title;
  alignment.tName = target.title;
  int insertion = 5;
  int deletion  = 5;
  DistanceMatrixScoreFunction<FASTASequence, FASTASequence> distScoreFn(SMRTDistanceMatrix, insertion, deletion);
  
  OneGapAlign(query, 
              leftTarget, 
              rightTarget, 
              distanceBetweenLeftAndRight,
              distScoreFn, alignment);

  StickPrintAlignment(alignment, query, target, cout);
}
