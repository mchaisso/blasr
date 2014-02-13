#ifndef POSITION_TABLE_H_
#define POSITION_TABLE_H_

#include <string>
using namespace std;

class PositionTable {
 public:
	int tableLength;
	int **table;
	int *rowLength;
	int InitFromFile(string &fileName);
};





#endif
