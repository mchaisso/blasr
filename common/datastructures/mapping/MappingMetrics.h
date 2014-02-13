#ifndef DATASTRUCTURES_MAPPING_MAPPING_METRICS_H_
#define DATASTRUCTURES_MAPPING_MAPPING_METRICS_H_

#include <iostream>
#include <time.h>
#include <map>
#ifdef __APPLE__
#pragma weak clock_gettime
#include <mach/mach.h>
#include <mach/clock.h>
#include <mach/mach_time.h>
#include <errno.h>
typedef enum {
  CLOCK_REALTIME,
  CLOCK_MONOTONIC,
  CLOCK_PROCESS_CPUTIME_ID,
  CLOCK_THREAD_CPUTIME_ID
} clockid_t;
static mach_timebase_info_data_t __clock_gettime_inf;

int clock_gettime(clockid_t clk_id, struct timespec *tp) {
  kern_return_t   ret;
  clock_serv_t    clk;
  clock_id_t clk_serv_id;
  mach_timespec_t tm;
  
  uint64_t start, end, delta, nano;
  
  task_basic_info_data_t tinfo;
  task_thread_times_info_data_t ttinfo;
  mach_msg_type_number_t tflag;
  
  int retval = -1;
  switch (clk_id) {
  case CLOCK_REALTIME:
  case CLOCK_MONOTONIC:
    clk_serv_id = clk_id == CLOCK_REALTIME ? CALENDAR_CLOCK : SYSTEM_CLOCK;
    if (KERN_SUCCESS == (ret = host_get_clock_service(mach_host_self(), clk_serv_id, &clk))) {
      if (KERN_SUCCESS == (ret = clock_get_time(clk, &tm))) {
	tp->tv_sec  = tm.tv_sec;
	tp->tv_nsec = tm.tv_nsec;
	retval = 0;
      }
    }
    if (KERN_SUCCESS != ret) {
      errno = EINVAL;
      retval = -1;
    }
    break;
  case CLOCK_PROCESS_CPUTIME_ID:
  case CLOCK_THREAD_CPUTIME_ID:
    start = mach_absolute_time();
    if (clk_id == CLOCK_PROCESS_CPUTIME_ID) {
      getpid();
    } else {
      sched_yield();
    }
    end = mach_absolute_time();
    delta = end - start;
    if (0 == __clock_gettime_inf.denom) {
      mach_timebase_info(&__clock_gettime_inf);
    }
    nano = delta * __clock_gettime_inf.numer / __clock_gettime_inf.denom;
    tp->tv_sec = nano * 1e-9;  
    tp->tv_nsec = nano - (tp->tv_sec * 1e9);
    retval = 0;
    break;
  default:
    errno = EINVAL;
    retval = -1;
  }
  return retval;
}
#endif // __APPLE__

class Timer {
 public:
	bool keepHistogram, keepList;
	timespec cpuclock[2];
	int elapsedClockMsec;
	float   elapsedTime;
	map<int,int> histogram;
  vector<int> msecList;
	long long totalElapsedClock;
  string header;
  
  
	Timer(string _header="") {
		keepHistogram = false;
    keepList      = false;
		totalElapsedClock = 0;
    header        = _header;
    elapsedClockMsec = 0;
    elapsedTime   = 0.0;
	}

  int ListSize() {
    return msecList.size();
  }

  void PrintHeader(ostream &out) {
    if (msecList.size() > 0) {
      out << header << " ";
    }
  }

  void PrintListValue(ostream &out, int index) {
    if (msecList.size() > 0) {
      out << msecList[index] << " ";
    }
  }
	void Tick() {
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &cpuclock[0]);
	}

  void SetStoreElapsedTime(bool value) {
    keepList = value;
  }

  void SetStoreHistgram(bool value) {
    keepHistogram = value;
  }
    
	void Tock() {
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &cpuclock[1]);
		elapsedClockMsec   = (cpuclock[1].tv_nsec - cpuclock[0].tv_nsec)/1000;
		totalElapsedClock += elapsedClockMsec;
		elapsedTime  = ((1.0)*elapsedClockMsec);
		if (keepHistogram) {
			// keep a histogram in number of milliseconds per operation
			if (histogram.find(elapsedClockMsec) == histogram.end()) {
				histogram[elapsedClockMsec] = 1;
			}
			else {
				histogram[elapsedClockMsec]++;
			}
		}
    if (keepList) {
      msecList.push_back(elapsedClockMsec);
    }
	}

	void Add(const Timer &rhs) {
		elapsedClockMsec += rhs.elapsedClockMsec;
		elapsedTime  += rhs.elapsedTime;
		totalElapsedClock += rhs.totalElapsedClock;
    msecList.insert(msecList.end(), rhs.msecList.begin(), rhs.msecList.end());
	}

  void SetHeader(string _header) {
    header = _header;
  }

};


class MappingClocks {
 public:
	Timer total;
	Timer findAnchors;
	Timer mapToGenome;
	Timer sortMatchPosList;
	Timer findMaxIncreasingInterval;
	Timer alignIntervals;           
  vector<int> nCellsPerSample;
  vector<int> nBasesPerSample;

  void AddCells(int nCells) {
    nCellsPerSample.push_back(nCells);
  }

  void AddBases(int nBases) {
    nBasesPerSample.push_back(nBases);
  }
    
  int  GetSize() {
    return total.ListSize();
  }

  MappingClocks() {
    total.SetHeader("Total");
    findAnchors.SetHeader("FindAnchors");
    mapToGenome.SetHeader("MapToGenome");
    sortMatchPosList.SetHeader("SortMatchPosList");
    findMaxIncreasingInterval.SetHeader("FindMaxIncreasingInterval");
    alignIntervals.SetHeader("AlignIntervals");
  }

  void PrintHeader(ostream &out) {
    total.PrintHeader(out);
    findAnchors.PrintHeader(out);
    mapToGenome.PrintHeader(out);
    sortMatchPosList.PrintHeader(out);
    findMaxIncreasingInterval.PrintHeader(out);
    alignIntervals.PrintHeader(out); 
  }

  void PrintList(ostream &out, int index) {
    
    total.PrintListValue(out,index);
    findAnchors.PrintListValue(out,index);
    mapToGenome.PrintListValue(out,index);
    sortMatchPosList.PrintListValue(out,index);
    findMaxIncreasingInterval.PrintListValue(out,index);
    alignIntervals.PrintListValue(out,index);
    if (nCellsPerSample.size() > 0) {
      out << nCellsPerSample[index] << " ";
    }
    if (nBasesPerSample.size() > 0) {
      out << nBasesPerSample[index] << " ";
    }
    out << endl;
  }

  void SetStoreList(bool value=true) {
    total.SetStoreElapsedTime(value);
    findAnchors.SetStoreElapsedTime(value);
    mapToGenome.SetStoreElapsedTime(value);
    sortMatchPosList.SetStoreElapsedTime(value);
    findMaxIncreasingInterval.SetStoreElapsedTime(value);
    alignIntervals.SetStoreElapsedTime(value);
  }

	void AddClockTime(const MappingClocks &rhs) {
		total.Add(rhs.total);
		findAnchors.Add(rhs.findAnchors);
		mapToGenome.Add(rhs.mapToGenome);
		sortMatchPosList.Add(rhs.sortMatchPosList);
		findMaxIncreasingInterval.Add(rhs.findMaxIncreasingInterval);
		alignIntervals.Add(rhs.alignIntervals);
	}
};


class MappingMetrics {
 public:
	MappingClocks clocks;
	int numReads;
	int numMappedReads;
  int numMappedBases;
  vector<int> mappedBases;
  vector<int> cellsPerAlignment;
  vector<int> anchorsPerAlignment;
  vector<int> sdpAnchors, sdpBases, sdpClock;
	long totalAnchors;
	int anchorsPerRead;
	long totalAnchorsForMappedReads;
	
	MappingMetrics() {
		numReads = 0;
		numMappedReads = 0;
    numMappedBases = 0;
    anchorsPerRead = 0;
		totalAnchorsForMappedReads = 0;
		totalAnchors = 0;
	}
  void StoreSDPPoint(int nBases, int nSDPAnchors, int nClock) {
    sdpBases.push_back(nBases);
    sdpAnchors.push_back(nSDPAnchors);
    sdpClock.push_back(nClock);
  }

  void SetStoreList(bool value=true) {
    clocks.SetStoreList(value);
  }

	void PrintSeconds(ostream&out, long sec ){
		out << sec << " Msec";
	}

	void PrintFraction(ostream &out, float frac) {
		out << setprecision(2) << frac;
	}

  void RecordNumAlignedBases(int nBases) {
    mappedBases.push_back(nBases);
  }

  void RecordNumCells(int nCells) {
    cellsPerAlignment.push_back(nCells);
  }

	void Collect(MappingMetrics &rhs) {
		clocks.AddClockTime(rhs.clocks);
		totalAnchors += rhs.totalAnchors;
		numReads += rhs.numReads;
		numMappedReads += rhs.numMappedReads;
		totalAnchorsForMappedReads += rhs.totalAnchorsForMappedReads;
    mappedBases.insert(mappedBases.end(), rhs.mappedBases.begin(), rhs.mappedBases.end());
    cellsPerAlignment.insert(cellsPerAlignment.end(), rhs.cellsPerAlignment.begin(), rhs.cellsPerAlignment.end());
	}

  void CollectSDPMetrics(MappingMetrics &rhs) {
    sdpAnchors.insert(sdpAnchors.end(), rhs.sdpAnchors.begin(), rhs.sdpAnchors.end());
    sdpBases.insert(sdpBases.end(), rhs.sdpBases.begin(), rhs.sdpBases.end());
    sdpClock.insert(sdpClock.end(), rhs.sdpClock.begin(), rhs.sdpClock.end());
  }

  void PrintSDPMetrics(ostream &out) {
    out << "nbases ncells time" << endl;
    int i;
    for (i = 0; i < sdpAnchors.size(); i++) {
      out << sdpBases[i] << " " << sdpAnchors[i] << " " << sdpClock[i] << endl;
    }
  }
    
  void PrintFullList(ostream &out) {
    // 
    // Print the full header
    //
    clocks.PrintHeader(out);
    out << " MappedBases Cells " << endl;
    //
    // Print all values: clocks + bases and cells.
    //
    int i;
    for (i = 0; i < clocks.GetSize(); i++) {
      clocks.PrintList(out,i);
      //      out << mappedBases[i] << " " << cellsPerAlignment[i] << endl;
    }
  }

	void PrintSummary(ostream &out) {
		out << "Examined " << numReads << endl;
		out << "Mapped   " << numMappedReads << endl;
		out << "Total mapping time\t";
		PrintSeconds(out, clocks.total.elapsedClockMsec);
		out << " \t";
		PrintSeconds(out, (1.0*clocks.total.elapsedClockMsec)/numReads);
		out << " /read" << endl;
		out << "      find anchors\t";
		PrintSeconds(out, clocks.mapToGenome.elapsedClockMsec);
		out << " \t";
		PrintSeconds(out, (1.0*clocks.mapToGenome.elapsedClockMsec)/numReads);
		out << endl;
		out << "      sort anchors\t";
		PrintSeconds(out, clocks.sortMatchPosList.elapsedClockMsec);
		out << " \t";
		PrintSeconds(out, (1.0*clocks.sortMatchPosList.elapsedClockMsec)/numReads);
		out << endl;
		out << " find max interval\t";
		PrintSeconds(out, clocks.findMaxIncreasingInterval.elapsedClockMsec);
		out << " \t";
		PrintSeconds(out, (1.0*clocks.findMaxIncreasingInterval.elapsedClockMsec)/numReads);
		out << endl;
		out << "Total anchors: " << totalAnchors << endl;
		out << "   Anchors per read: " << (1.0*totalAnchors) / numReads << endl;
		out << "Total mapped: " << totalAnchorsForMappedReads << endl;
		out << "   Anchors per mapped read: " << (1.0*totalAnchorsForMappedReads) / numMappedReads << endl;
	}
	
	void AddClock(MappingClocks &clocks) {
		clocks.AddClockTime(clocks);
	}
};



#endif
