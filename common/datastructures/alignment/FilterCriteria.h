/*
 * =====================================================================================
 *
 *       Filename:  FilterCriteria.h
 *
 *    Description:  Criteria for filtering alignment hits. 
 *
 *        Version:  1.0
 *        Created:  03/22/2013 11:33:00 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#ifndef FILTER_CRITERIA_H
#define FILTER_CRITERIA_H
#include "datastructures/alignment/AlignmentCandidate.h"

enum SCORESIGN {NEG=-1, POS=1};

class Score {
public:
    double    score;
    SCORESIGN sign;

    Score() {score = 0; sign = NEG;}

    Score(const double & d, SCORESIGN ss) {
        score = d;
        sign  = ss;
    }

    Score(const double & d, int & ss) {
        score = d;
        SetScoreSign(ss);
    }

    SCORESIGN SetScoreSign(int & ss) {
        if (ss == POS) 
            sign = POS;
        else if (ss == NEG) 
            sign = NEG;
        else {
            cout << "ERROR, invalid score sign " << ss << endl;
            exit(1);
        }
        return sign;
    }

    // True if two scores are equivalent.
    bool Equal(const double & d) {
        double delta = d - score;
        double errorunit = 1e-10;
        if (delta < errorunit && delta > -errorunit) 
            return true;
        return false;
    }

    // POS, the greater the better.
    // NEG, the less the better.
    bool BetterThan(const double & d) {
        if (Equal(d)) return false;
        if (sign == POS) return (score > d);
        else return (score < d);
    }

    bool BetterThanOrEqual(const double & d) {
        return BetterThan(d) or Equal(d);
    }

    bool WorseThan(const double & d) {
        return not (BetterThanOrEqual(d));
    }



};


class FilterCriteria {
public:
    int minAccuracy;
    int minLength;
    int minPctSimilarity;
    Score scoreCutoff; 
    bool useScore;
    int verbosity;

    // Potential fileters to use 
    int minAnchorBases;
    int minAnchorSize;
    int minZ;

    FilterCriteria() {
        SetDefault();
    }

    // Set default filter criteria.
    void SetDefault() {
        minAccuracy   = 70;
        minLength = 50;
        minPctSimilarity = 70;
        minAnchorSize = 12;
        minZ          = 0;
        useScore      = false; 
        verbosity     = 0;
        // score is not a filter criteria unless sepecified.
    }

    // Set score value.
    void SetScoreCutoff(double d) {
        scoreCutoff.score = d;
        useScore = true;
    }

    // Set score sign and return the sign.
    SCORESIGN SetScoreSign(int & ss) {
        return scoreCutoff.SetScoreSign(ss);
    }

    SCORESIGN SetScoreSign(SCORESIGN ss) {
        scoreCutoff.sign = ss;
        return ss;
    }

    void SetMinAccuracy(const int & a) {
        minAccuracy = a;
    }
    void SetMinReadLength(const int & l) {
        minLength = l;
    }
    void SetMinPctSimilarity(const int & s) {
        minPctSimilarity = s;
    }
    void SetMinAnchorSize(const int & a) {
        minAnchorSize = a;
    }
    void SetMinZ(const int & z) {
        minZ = z;
    }
    
    // Return whether or not an alignment satisfies all filter criteria
    bool Satisfy(AlignmentCandidate<> & alignment) {
    //
    // tAlignedSeq, qAlignedSeq, tAlignedSeqPos, qAlignedSeqPos
    // tAlignedSeqLength, qALignedSeqLength, readIndex, tIndex, 
    // mapQV, clusterScore, clusterWeight, tIsSbustring, qIsSubstring
    // tTitle, qTitle, nMatch, nIns, nDel, pctSimilarity
    // 
        if (alignment.qAlignedSeqLength   < minLength) {
            if (verbosity > 0) 
                cout << "Alignment length is too short (" 
                     << alignment.qAlignedSeqLength 
                     << " < " << minLength << ")." <<endl;
            return false;
        } 
        if (alignment.pctSimilarity * 100 < minPctSimilarity) {
            if (verbosity > 0) 
                cout << "Percentage similarity (" 
                     << alignment.pctSimilarity
                     << " < " << minPctSimilarity 
                     << ") is too low." << endl;
            return false;
        }
        
        if (useScore && scoreCutoff.BetterThanOrEqual(alignment.score)) {
            if (verbosity > 0)
                cout << "Alignment score (" << alignment.score 
                     << " worse than " << scoreCutoff.score 
                     << ") is bad." << endl;
            return false;
        }

        double accuracy =  100 * alignment.nMatch / (double)(
                           alignment.nMismatch + alignment.nMatch + 
                           alignment.nIns + alignment.nDel);

        if (accuracy < minAccuracy) {
            if (verbosity > 0)
                cout << "Accuracy (" << accuracy << " < "
                     << minAccuracy << ") is too low." << endl;
            return false;
        }


        //int nAnchors, nAnchorBases;
        // nAnchors = nAnchorBases = 0;
        //Update alignment anchor info
        //alignment.ComputeNumAnchors(minAnchorSize, nAnchors, nAnchorBases);
        //if (nAnchors <= 0) {
        //    if (verbosity > 0) 
        //        cout << "No anchor exists in the alignment." << endl;
        //    return false;
        // }

        return true;
    }
};

#endif
