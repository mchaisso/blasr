#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include "NucConversion.h"

int main(int argc, char* argv[]) {

	ifstream bwtIn;
	bwtIn.open(argv[1]);
	/*
	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, bwtFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, bwtFile);
	bwtLength = BWTFileSizeInWord(bwt->textLength);
	fwrite(bwt->bwtCode, sizeof(unsigned int), bwtLength, bwtFile);
	fclose(bwtFile);
	*/
	unsigned int inverseSa0;
	unsigned int cumulativeFreq[4];
	unsigned int bwtLength;
	unsigned int *bwtCode;
	bwtIn.read((char*) &inverseSa0, sizeof(unsigned int));
	bwtIn.read((char*) cumulativeFreq, sizeof(unsigned int)*4);
	bwtIn.read((char*) &bwtLength, sizeof(unsigned int));
	unsigned int bwtWordLength = 	 (bwtLength + 16 - 1) / 16;
  bwtCode = new unsigned int[bwtWordLength];
	bwtIn.read((char*) bwtCode, sizeof(unsigned int ) * bwtWordLength);
	unsigned int w;
	unsigned int i;
	unsigned int bwtWord;
	w = 0;
	for (i = 0; i < bwtLength; i++ ){
		if (i % 16 == 0) {
			bwtWord = bwtCode[w++];
		}
		cout << TwoBitToAscii[bwtWord & 3];
		if (i % 50 == 49) 
			cout << endl;
		bwtWord = bwtWord >> 2;
	}
}
