#include "PositionTable.h"


#include <iostream>
#include <fstream>

int PositionTable::InitFromFile(string &fileName) {
	ifstream tableIn;
	tableIn.open(fileName.c_str(), ios_base::binary);
	if (!tableIn.good())
		return 0;
	
	tableIn.read((char*)&tableLength, sizeof(int));
	/* 
	  Init ragged arrays.
	*/

	if (tableLength == 0 ) {
		table = NULL;
		rowLength = NULL;
	}
	
	table = new int*[tableLength];
	rowLength = new int[tableLength];	

	//
	// Read in the array as a ragged array.
  //
	int i;
	for (i = 0; i < tableLength; i++) {
		tableIn.read((char*)&(rowLength[i]), sizeof(int));
		table[i] = new int[rowLength[i]];
		tableIn.read((char*)table[i], sizeof(int) * rowLength[i]);
	}
	return 1;
}
