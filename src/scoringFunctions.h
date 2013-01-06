#ifndef SCORING_FUNCTIONS_H
#define SCORING_FUNCTIONS_H

#include <vector>
#include "dp.h"

double gapPenalty(int fragSize, const AlignmentParams& ap);

double contigMissedSitePenalty(int dToClosestSite, const AlignmentParams& ap);

double scoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag,
                       const AlignmentParams& ap);

double scoringFunction2( const vector<FragData>::const_iterator& cB,
                         const vector<FragData>::const_iterator& cE,
                         const vector<FragData>::const_iterator& oB,
                         const vector<FragData>::const_iterator& oE,
                         bool boundaryFrag ,
                         const AlignmentParams& ap);

double localScoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag, const AlignmentParams& ap);


#endif
