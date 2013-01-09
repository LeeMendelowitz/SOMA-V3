#include "scoringFunctions.h"
#include "AlignmentParams.h"
#include "globals.h"

#include <cassert>

using namespace std;

// Penalty for "losing" a small contig fragment.
// (gapped alignment)
double gapPenalty(int fragSize, const AlignmentParams& ap)
{
    if (fragSize < ap.smallFrag)
        return -fragSize * ap.smallFragSlope;
    else
        return -Constants::INF;
}

// Penalty for missing a contig restriction site
double contigMissedSitePenalty(int dToClosestSite, const AlignmentParams& ap)
{
    if (dToClosestSite > ap.smallFrag)
        return -ap.C_r_contig;
    else
    {   
        // Add a small amount to the penalty to break ties with the gapPenalty
        // in boundary fragments where there is no sizing error
        return (gapPenalty(dToClosestSite, ap) + 1.0E-6);
    }
}


// nContigSites: number of unaligned contig sites
// nOpticalSites: number of unaligned optical sites
double scoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag,
                       const AlignmentParams& ap)
{
    // Model assumes that Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size 

    double chi2;
    
    if (boundaryFrag && contigLength < opticalLength)
        chi2 = 0.0; // Do not penalize boundary fragments for being too small
    else
    {
        double var = contigLength*ap.sigma2;
        double dl = (contigLength-opticalLength);
        chi2 = dl*dl/var;
        assert (chi2 >= -1E-12);
        // Check if this alignment block exceeds the maximum allowed sizing error
        if (chi2 > ap.chi2Max)
            return -Constants::INF;
    }
    return -(nContigSites*ap.C_r_contig +
           nOpticalSites*ap.C_r_optical +
           chi2);
}


// Scoring function which accounts for the distance to the closest
// aligned site for each unaligned interior site.
double scoringFunction2( const vector<FragData>::const_iterator& cB,
                         const vector<FragData>::const_iterator& cE,
                         const vector<FragData>::const_iterator& oB,
                         const vector<FragData>::const_iterator& oE,
                         bool boundaryFrag ,
                         const AlignmentParams& ap)
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
        double var = contigLength*ap.sigma2;
        double dl = (contigLength-opticalLength);
        chi2 = dl*dl/var;
        assert (chi2 >= -1E-12);
        // Check if this alignment block exceeds the maximum allowed sizing error
        if (chi2 > ap.chi2Max)
            return -Constants::INF;
    }
    chi2 = -1.0*chi2;

    // Compute miss penalty
    double contigMissScore = 0.0;
    int numOpticalFrags = oE - oB;
    double opticalMissScore = -ap.C_r_optical*(numOpticalFrags-1);
    int pos = 0;
    for(ci = cB; ci != cE-1; ci++)
    {
        pos += ci->size_;
        // The miss score is a function of the distance to the first missed site
        contigMissScore += contigMissedSitePenalty(min(pos, contigLength-pos), ap);
    }

    return (contigMissScore + opticalMissScore + chi2);
}

// Scoring functions for local alignment
// The greater (i.e. more positive) the score, the better the alignment
// nContigSites: number of unaligned contig sites
// nOpticalSites: number of unaligned optical sites
// contigLength: length of the contig alignment chunk
// opticalLength: length of the optical alignment chunk
// Returns a positive number for good alignments, a negative
// number for bad alignments

// The sizing score is given by a quadratic:
// f(delta) = (H/T^2)*(delta^2 - T^2)
// where delta = (opticalLength - contigLength)/VAR[opticalLength]
// Assuming that the observed optical fragment corresponds to the contig fragment,
// VAR[opticalLength] depends on the length of the contig chunk.
// Here, we use the model: VAR[opticalLength] = sigma^2*contigLength

double localScoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag, const AlignmentParams& ap)
{
    // Model assumes that Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size 
    double delta = (opticalLength - contigLength);
    double var = contigLength*ap.sigma2;
    double chi2 = delta*delta/var;

    // Check if this alignment block exceeds the maximum allowed sizing error
    if (chi2 > ap.chi2Max)
        return -Constants::INF;

    double sizeScore = -ap.A*(chi2 - ap.T2); // This is positive for small delta
    double missPenalty = -nContigSites*ap.C_r_contig - nOpticalSites*ap.C_r_optical;
    return sizeScore + missPenalty;
}
