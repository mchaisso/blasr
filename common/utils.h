#ifndef UTILS_H_
#define UTILS_H_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

/*
template <typename T>
int min(const T &a, const T&b) {
	return a < b;
}
*/



using namespace std;
template<typename t_file>
void CrucialOpen(string &fileName, t_file &file, std::ios_base::openmode mode= (std::ios_base::openmode)0 ) {
	if (mode==0)
		file.open(fileName.c_str());
	else
		file.open(fileName.c_str(), mode);

	if (!file.good()) {
		cout << "Could not open " << fileName << endl;
		exit(1);
	}
}
template<typename T_Int>
T_Int CeilOfFraction(T_Int num, T_Int denom) {
	return num / denom + ((num % denom) && 1);
}

#endif
