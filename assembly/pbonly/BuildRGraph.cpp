#include "RGraph.h"
#include "ReadNode.h"
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cout << "usage: testRGraph file.rm4 " << endl;
    exit(1);
  }
  string rm4FileName = argv[1];
  RGraph rGraph;
  rGraph.BuildFromRM4(rm4FileName);

  //  rGraph.PrintPossibleChimeras(2);
  //  rGraph.PrintLowCoverageReads(1000000);

  cout << "There are " << ReadNode::counter << " nodes and " << ReadNodeSuperset::counter << " supernodes." << endl;
  rGraph.RemoveLowCoverage(5);
  cout << "Removed " <<         ReadNodeSuperset::countDeleted << " out of " << ReadNodeSuperset::counter << endl;
  
  rGraph.MergeOverlaps();

  int nConsistent = rGraph.CountConsistentSupernodes();
  cout << nConsistent << " out of " <<  ReadNodeSuperset::counter - ReadNodeSuperset::countDeleted << endl;
  
  return 0;
}
