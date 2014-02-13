/*
 * ============================================================================
 *
 *       Filename:  RangeUtils.h
 *
 *    Description:  Parse a list of ranges separated by commas.
 *                  Designed for specifying a set of holeNumbers to analyze.
 *
 *        Version:  1.0
 *        Created:  05/02/2013 12:51:27 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ============================================================================
 */
#include <stdlib.h>
#include <algorithm>
#include "StringUtils.h"
#include "../Types.h"

class Range {
public:
    UInt start;  // Start position of a range, inclusive
    UInt end;    // End position of a range, inclusive
    Range(UInt pStart) {
        start = end = pStart;
    }
    Range(UInt pStart, UInt pEnd) {
        start = pStart;
        end = pEnd;
        if (start > end) {
            cout << "ERROR: start of a range should be less than the end."
                 << endl;
            exit(1);
        }
    }
    bool contains(const UInt & query) {
        return (start <= query && query <= end);
    }
    bool operator < (const Range & pRange) const {
        if (start == pRange.start) {
            return (end > pRange.end);
        }
        return (start < pRange.start);
    }
};

// Input is a comma-delimited string of ranges.
// e.g. 1,2,3,10-20
bool ParseRanges(string & rangesStr, 
                 vector<Range> & ranges) {
    ranges.clear();
    bool parseSucceed = true;
    vector<string> strList;
    ParseSeparatedList(rangesStr, strList, ',');
    for(int i=0; i<int(strList.size()); i++) {
        string & str = strList[i];
        if(str.find('-') == string::npos) {
            ranges.push_back(Range(atoi(str.c_str())));
        } else {
            vector<string> start_end;
            ParseSeparatedList(str, start_end, '-');
            if (start_end.size() != 2) {
                parseSucceed = false;
                break;
            }
            ranges.push_back(Range(atoi(start_end[0].c_str()),
                        atoi(start_end[1].c_str())));
        }
    }
    if (parseSucceed) {
        sort(ranges.begin(), ranges.end());
    } else {
        ranges.clear();
    }
    return parseSucceed;
}

class Ranges {
public:
    vector<Range> ranges;
    Ranges(string rangesStr="") {
        setRanges(rangesStr);
    }
    bool setRanges(string rangesStr) {
        return ParseRanges(rangesStr, ranges);
    }
    
    int size() {
        return ranges.size();
    }
    UInt max() {
        if (size() == 0) {
            cout << "ERROR, could not determine the maximum value "
                 << "of an empty Ranges object." << endl;
            exit(1);
        }
        return ranges.back().end;
    }

    bool contains(const UInt & query) {
        vector<Range> searchRanges;
        searchRanges.push_back(Range(0, ranges.size()-1));
        while (searchRanges.size() > 0) {
            Range searchRange = searchRanges.back();
            searchRanges.pop_back();
            UInt mid = (searchRange.start + searchRange.end) / 2;
            if (ranges[mid].contains(query)) {
                return true;
            }
            if (mid > 0 && searchRange.start <= mid - 1) {
                searchRanges.push_back(Range(searchRange.start,
                                       mid - 1));
            }
            if (ranges[mid].start <= query and
                searchRange.end >= mid + 1) {
                searchRanges.push_back(Range(mid + 1,
                                       searchRange.end));
            }
        }
        return false;
    }
};


