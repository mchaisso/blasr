#include "data/hdf/HDFBasReader.h"
#include "SMRTSequence.h"
#include <string>

int main(int argc, char* argv[]) {
	
	string basFileName = argv[1];

  HDFBasReader reader;
  reader.InitializeDefaultIncludedFields();
	reader.Initialize(basFileName);

	SMRTSequence read;
  BaseFile baseFile;
  reader.ReadBaseFile(baseFile);
	return 0;
  
}


