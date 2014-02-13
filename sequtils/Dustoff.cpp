#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include "../common/utils.h"
using namespace std;
int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "usage: dustoff sequence" << endl;
    cout << "Prints the locations of non 'ACGT' nucleotides." << endl;
		exit(1);
	}
	string inFile = argv[1];
	ifstream in;
	ofstream out;

	CrucialOpen(inFile, in);
	int lineNumber = 0;
	char ch;
	int linePos = 0;
  cout << "line_number line_pos (int) char" << endl;
	while (in) {
		ch = toupper(in.get());
		if (ch == '>') {
			string line;
			getline(in,line);
      ++lineNumber;
      linePos = 0;
		}
    else {
      if (ch != 'A' and ch != 'C' and ch != 'G' and ch != 'T' and ch != 'N' and ch != '\n' and ch != '\0') {
        cout << lineNumber << " " << linePos << " " << (int) ch << " " << ch << endl;
      }
      if (ch == '\n') {
        ++lineNumber;
        linePos = 0;
      }
      else {
        linePos++;
      }
    }
	}
}

	
