#ifndef GLOBAL_SCORER_H
#define GLOBAL_SCORER_H

#include <cassert>

#include "Scorer.h"
#include "globals.h"

// In this scoring scheme, the scores are negative (representing penalties).
// The more positive (less negative) the score is, the better.

class GlobalScorer
{

    public:
    GlobalScorer(const AlignmentParams& ap) :
        ap_(ap) {}
    ~GlobalScorer() {}
    inline Score scoreGap(int gapSize) const ;
    inline Score scoreAlignment(const std::vector<FragData>::const_iterator cB,
                         const std::vector<FragData>::const_iterator cE,
                         const std::vector<FragData>::const_iterator oB,
                         const std::vector<FragData>::const_iterator oE,
                         bool boundaryFrag) const;
    void scoreMatchResult(MatchResult * pResult) const;
    Score scoreMatchedChunk(MatchedChunk& chunk) const;

    private:
    double contigMissedSitePenalty(int dToClosestSite) const;
    AlignmentParams ap_;
};

// Penalty for "losing" a small contig fragment.
// (gapped alignment)
inline Score GlobalScorer::scoreGap(int gapSize) const
{
        double contigScore;
        if (gapSize < ap_.smallFrag)
            contigScore =  -gapSize * ap_.smallFragSlope;
        else
            contigScore = -Constants::INF;
        return Score(contigScore, 0.0, 0.0);
}

// Penalty for missing a contig restriction site
inline double GlobalScorer::contigMissedSitePenalty(int dToClosestSite) const
{
    if (dToClosestSite > ap_.smallFrag)
    {
        return -ap_.C_r_contig;
    }
    else
    {   
        // Add a small amount to the penalty to break ties with the gapPenalty
        // in boundary fragments where there is no sizing error
        Score score = scoreGap(dToClosestSite);
        return score.contig += 1.0E-6;
    }
}


// Scoring function which accounts for the distance to the closest
// aligned site for each unaligned interior site.
inline Score GlobalScorer::scoreAlignment( const vector<FragData>::const_iterator cB,
                         const vector<FragData>::const_iterator cE,
                         const vector<FragData>::const_iterator oB,
                         const vector<FragData>::const_iterator oE,
                         bool boundaryFrag) const
{
    vector<FragData>::const_iterator ci, oi;

    // Calculate total contig and optical block length
    int contigLength = 0;
    for (ci = cB; ci != cE; ci++)
        contigLength += ci->size_;

    int opticalLength = 0;
    for (oi = oB; oi != oE; oi++)
        opticalLength += oi->size_;

    // Compute sizing error
    double chi2 = 0.0;
    if (boundaryFrag && contigLength < opticalLength)
        chi2 = 0.0; // Do not penalize boundary fragments for being too small
    else
    {
        double var = contigLength*ap_.sigma2;
        double dl = (contigLength-opticalLength);
        chi2 = dl*dl/var;
        assert (chi2 >= -1E-12);
        // Check if this alignment block exceeds the maximum allowed sizing error
        if (chi2 > ap_.chi2Max)
            return Score(0.0, 0.0, -Constants::INF);
    }
    chi2 = -1.0*chi2;

    // Compute miss penalty
    double contigMissScore = 0.0;
    int numOpticalFrags = oE - oB;
    double opticalMissScore = -ap_.C_r_optical*(numOpticalFrags-1);
    int pos = 0;
    for(ci = cB; ci != cE-1; ci++)
    {
        pos += ci->size_;
        // The miss score is a function of the distance to the first missed site
        contigMissScore += contigMissedSitePenalty(min(pos, contigLength-pos));
    }

    return Score(contigMissScore, opticalMissScore, chi2);
}

#endif
