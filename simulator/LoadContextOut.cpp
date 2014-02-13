#include "../common/simulator/ContextOutputList.h"
#include "../common/utils.h"
#include <string>
using namespace std;


int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "usage: loadContextOut in" << endl;
		exit(1);
	}
	string inFileName;
	inFileName = argv[1];
	
	ifstream tableIn;
	CrucialOpen(inFileName, tableIn, std::ios::in);

	ContextOutputList outputList;
	outputList.Read(tableIn);
	int nContextElements = 1 << (2*outputList.contextLength);
	cout << outputList.outputMap.size() << " contexts stored of " << nContextElements << " " << endl;
	return 0;
}
