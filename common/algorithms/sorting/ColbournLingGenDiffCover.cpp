#include <iostream>
#include <stdlib.h>
#include <vector>
using namespace std;
int main(int argc, char* argv[]) {
	
	int r = atoi(argv[1]);

	int d = 6*r+4;

	vector<int> B;
	B.resize(d-1);
	
	int i;
	int bi = 0;
	cout << d << endl;
	// Fill from the terms:
	// 1^d(d+1)^1(2*d+1)^d(4*d+3)^{2*r+1}(2*d+2)^{d+1}1^d
	// term 1
	for (i = 0; i < r; i++) {
		B[bi] = 1;
		bi++;
	}
	// term 2
	B[bi] = (r+1);
	bi++;
	for (i = 0; i < r; i++) {
		B[bi] = 2*r+1;
		bi++;
	}

	// term 3
	for (i = 0; i < 2*r+1; i++) {
		B[bi] = 4*r+3;
		bi++;
	}

	// term 4:
	for (i = 0; i < r+1; i++) {
		B[bi] = 2*r+2;
		bi++;
	}
	
	// term 5: 1
	for (i = 0; i < r; i++ ){
		B[bi] = 1;
		bi++;
	}
	int a = 0;
	cout << 24*r*r+36*r+13 << endl;
	for (i = 0; i < d; i++ ) {
		cout << a << ",";
		a = a + B[i];
	}
}
	
