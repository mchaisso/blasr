#include "../common/datastructures/metagenome/TitleTable.h"
#include "../common/FASTAReader.h"
#include <string>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
  
  if (argc < 3) {
    cout << "usage: printTitleTable fastaFile titleTable" << endl;
    exit(1);
  }
  string sourceFastaName = argv[1];
  string titleTableName  = argv[2];

  FASTAReader fastaReader;
  vector<string> titles;
  fastaReader.Init(sourceFastaName);
  
  ofstream out;
  CrucialOpen(titleTableName, out, std::ios::out);
  FASTASequence seq;
  while(fastaReader.GetNext(seq)) {
    titles.push_back(seq.title);
  }
  TitleTable titleTable;
  titleTable.CopyFromVector(titles);

  titleTable.Write(out);

  return 0;
}
