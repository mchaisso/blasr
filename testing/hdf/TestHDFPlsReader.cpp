#include "data/hdf/HDFPlsReader.h"
#include "SMRTSequence.h"
#include <string>

int main(int argc, char* argv[]) {
	
	string plsFile = argv[1];

  HDFPlsReader reader;
  reader.IncludeField("MidSignal");
  reader.IncludeField("MeanSignal");
  reader.IncludeField("MaxSignal");
	reader.Initialize(plsFile);

	SMRTSequence read;
  PulseFile pulseFile;
  reader.ReadPulseFile(pulseFile);
	return 0;
  
}


