#include "cdna/GencodeGFFFile.h"
#include "cdna/GencodeGFFGene.h"
#include "cdna/GeneDB.h"

#include <map>
#include <set>

using namespace std;

int main(int argc, char* argv[]) {
  string alignmentsFileName, gencodeDBFileName;
  
  if (argc < 2) {
    cout << "usage: validateExons alignmentsFile gencodeDBFile" << endl;
    exit(1);
  }

  alignmentsFileName = argv[1];
  gencodeDBFileName  = argv[2];

  GeneDB genedb;
  genedb.Read(gencodeDBFileName);
  

}
