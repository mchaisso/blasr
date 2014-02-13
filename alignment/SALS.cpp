#include "../common/datastructures/suffixarray/SuffixArray.h"
#include "../common/datastructures/suffixarray/SuffixArrayTypes.h"
#include "../common/utils.h"

using namespace std;
#include <string>
#include <iostream>
int main(int argc, char* argv[]) {

  if (argc <= 1) {
    cout << "usage sals genome.sa" <<endl;
    exit(1);
  }

  string saFileName = argv[1];
  
  DNASuffixArray sa;
  if (!sa.LightRead(saFileName)) {
    cout << "The file is not in a sa format." << endl;
    exit(1);
  }
  
  if (sa.componentList[DNASuffixArray::CompArray]) {
    cout << " * has a suffix array." << endl;
  }
  else {
    cout << " * does not contain a suffix array." << endl;
  }

  if (sa.componentList[DNASuffixArray::CompLookupTable]) {
    cout << " * has a lookup table for word size. " << sa.lookupPrefixLength 
         << endl;
  }  
  else {
    cout << " * does not have a lookup table." << endl;
  }
}

