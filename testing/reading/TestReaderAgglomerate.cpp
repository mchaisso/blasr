#include "../../common/files/ReaderAgglomerate.h"
#include "FASTASequence.h"
#include <string>


int main(int argc, char* argv[]) {
	
	string fileName = argv[1];
	int stride = 0;
	if (argc > 2) {
		stride = atoi(argv[2]);
	}

	
}
