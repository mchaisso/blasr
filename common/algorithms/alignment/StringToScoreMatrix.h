#ifndef STRING_TO_SCORE_MATRIX_H_
#define STRING_TO_SCORE_MATRIX_H_

#include <string>
#include <sstream>
#include <vector>
bool  StringToScoreMatrix(string &str, int matrix[5][5]) {
  stringstream strm(str);
  vector<int> values;
  while(strm) {
    int val;
    if ((strm >> val)) {
      values.push_back(val);
    }
  }
  if (values.size() != 25) {
    return 0;
  }
  else {
    int i,j;
    int index = 0;
    for (i = 0; i < 5; i++) {
      for (j = 0; j < 5; j++) {
        matrix[i][j] = values[index];
        ++index;
      }
    }
    return true;
  }
}


#endif
