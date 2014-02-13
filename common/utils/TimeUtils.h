#ifndef UTILS_TIME_UTILS_H_
#define UTILS_TIME_UTILS_H_

#include <ctime>
#include <sstream>
#include <time.h>
string GetTimestamp() {
  time_t timer;
  time(&timer);  // t is an integer type
  // Prepare timestamp in the format : 2012-04-05T09:26:02.689093
  stringstream timeStrm;
  struct tm t;
  localtime_r(&timer, &t);
  timeStrm << t.tm_year + 1900 << "-"
           << setfill('0') << setw(2) << t.tm_mon + 1 << "-"
           << setfill('0') << setw(2) << t.tm_mday << "T"
           << setfill('0') << setw(2) << t.tm_hour << ":"
           << setfill('0') << setw(2) << t.tm_min << ":"
           << setfill('0') << setw(2) << t.tm_sec;
  return timeStrm.str();
}

#endif
