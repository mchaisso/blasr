#include "RGraph.h"
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
}
