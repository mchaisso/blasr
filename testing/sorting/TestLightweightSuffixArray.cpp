#include "algorithms/sorting/LightweightSuffixArray.h"
#include <iostream>
#include <string>
using namespace std;
#include "Types.h"
int main(int argc, char* argv[]) {
	UInt D[] = {1,2,4};
	UInt lengthD = 3;
	UInt v = 7;

	char text[] = "a rose is a rose is a rose       ";
	UInt textLength = 26;
	UInt *index = new UInt[textLength+1];

	LightweightSuffixSort((unsigned char*) text, textLength, index, v);
	UInt i;
	for (i = 0; i < textLength; i++ ) {
		cout << i << " " << index[i] << " " << &text[index[i]] << endl;
	}

	return 0;
}
	
