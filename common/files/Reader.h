#ifndef READER_H_
#define READER_H_

class Reader {
 public:
  virtual bool HasRegionTable() {
    return false;
  }
  
  virtual void Initialize() {
    cout << " This should be over ridden! " << endl;
  }

  virtual void Advance(int nSteps) {
  }

};


#endif
