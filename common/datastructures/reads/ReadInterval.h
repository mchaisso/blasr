#ifndef DATASTRUCTURES_READS_READ_INTERVAL
#define DATASTRUCTURES_READS_READ_INTERVAL

class ReadInterval {
public:
	int start;
	int end;
  int score;
 ReadInterval(int s, int e, int sc=0) : start(s), end(e), score(sc) {};
	ReadInterval& operator=(const ReadInterval &rhs) {
		start = rhs.start;
		end   = rhs.end;
    score = rhs.score;
		return *this;
	}
};

#endif
