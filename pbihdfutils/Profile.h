#ifndef PBIHDFUTILS_PROFILE_H_
#define PBIHDFUTILS_PROFILE_H_

class Profile {
public:
  Matrix<int> profileMatrix;
  static const int delIndex = 4;
  static const int insOffset = 5;

  Profile(int length) {
    profileMatrix.Resize(9,length);
    profileMatrix.Initialize(0);
  }

  int IncrementNucCount(char nucleotide, int pos) {
    if (nucleotide == ' ') {
      profileMatrix[delIndex][pos]++;
    }
    else {
      if (nucleotide != '.') {
        assert(pos < profileMatrix.GetNCols());
        assert(TwoBit[nucleotide] < 4);
        profileMatrix[TwoBit[nucleotide]][pos]++;
      }
      else {
        cout << "ERROR! Cannot produce a profile when miscallView is used." << endl;
        exit(1);
      }
    }
  }

  int IncrementInsCount(char nucleotide, int pos) {
    assert(nucleotide != ' ');
    profileMatrix[TwoBit[nucleotide]+insOffset][pos]++;
  }

  void StoreProfile(string alignedString, int refStart, int readStart, vector<InsertedString> &insertions, int offset=0) {
    int i;
    int start, end;
    for (i = 0; i < alignedString.size() and alignedString[i] == ' '; i++);
    start = i;
    for (i = alignedString.size(); i > 0 and alignedString[i-1] == ' '; i--);
    end   = i;
    for (i = start; i < end; i++) {
      IncrementNucCount(alignedString[i], i+offset);
    }
    for (i = 0; i < insertions.size(); i++) {
      IncrementInsCount(insertions[i].insSeq[0], insertions[i].pos + readStart - refStart + offset);
    }
  } 

  void PrintRow(int rowIndex) {
    int i;
    if (profileMatrix.GetNCols() == 0) {
      return;
    }
    for (i = 0; i + 1 < profileMatrix.GetNCols(); i++) {
      cout.width(3);
      cout << profileMatrix[rowIndex][i]<< ",";
    }
    cout.width(3);
    cout << profileMatrix[rowIndex][i] << endl;
  }

  void Print(int marginWidth) {
    string padding;
    if (marginWidth < 2) {
      cout << "margin width must be greater than 2" << endl;
      assert(0);
    }
    padding.resize(marginWidth - 2);
    fill(padding.begin(), padding.end(), ' ');
    cout << " A" << padding;
    PrintRow(0);
    cout << " C" << padding;
    PrintRow(1);
    cout << " G" << padding;
    PrintRow(2);
    cout << " T" << padding;
    PrintRow(3);
    cout << " D" << padding;
    PrintRow(4);
    cout << "IA" << padding;
    PrintRow(5);
    cout << "IC" << padding;
    PrintRow(6);
    cout << "IG" << padding;
    PrintRow(7);
    cout << "IT" << padding;
    PrintRow(8);
  }
};

#endif
