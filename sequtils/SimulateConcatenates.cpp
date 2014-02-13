#include "FASTAReader.h"
#include "FASTASequence.h"
#include "CommandLineParser.h"
#include "statistics/statutils.h"
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
  
  CommandLineParser clp;

  clp.RegisterStringOption("");

  InitializeRandomGeneratorWithTime();

  

}
