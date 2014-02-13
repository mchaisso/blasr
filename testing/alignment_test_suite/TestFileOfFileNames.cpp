#include "utils/FileOfFileNames.h"
#include "CommandLineParser.h"

int main(int argc, char* argv[]) {

  commandLineParser clp;
  vector<string> fileNames;
  clp.RegisterStringListOption("files", &fileNames, "a list of file names", true);

  clp.ParseCommandLine(argc, argv);

  FileOfFileNames::FlattenFileNameList(fileNames);

  int f;
  for (f = 0; f < fileNames.size(); f++) {
    cout << fileNames[f] << endl;
  }
};
  
