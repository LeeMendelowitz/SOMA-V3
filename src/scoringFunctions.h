#ifndef SCORING_FUNCTIONS_H
#define SCORING_FUNCTIONS_H

#include <vector>
#include <cassert>
#include "mapTypes.h"
#include "globals.h"
#include "AlignmentParams.h"

double gapPenalty(int fragSize, const AlignmentParams& ap);

double contigMissedSitePenalty(int dToClosestSite, const AlignmentParams& ap);

double scoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag,
                       const AlignmentParams& ap);

double scoringFunction2( const std::vector<FragData>::const_iterator& cB,
                         const std::vector<FragData>::const_iterator& cE,
                         const std::vector<FragData>::const_iterator& oB,
                         const std::vector<FragData>::const_iterator& oE,
                         bool boundaryFrag ,
                         const AlignmentParams& ap);

double localScoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag, const AlignmentParams& ap);


// nContigSites: number of unaligned contig sites
// nOpticalSites: number of unaligned optical sites
inline double scoringFunction(int nContigSites, int nOpticalSites,
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

#endif
