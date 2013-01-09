#ifndef SCORING_FUNCTIONS_H
#define SCORING_FUNCTIONS_H

#include <vector>
#include "mapTypes.h"

//Forward declaration
class AlignmentParams;

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


#endif
