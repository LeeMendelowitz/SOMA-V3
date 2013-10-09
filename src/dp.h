// dp.h
#ifndef DP_H
#define DP_H

#include <vector>

class MatchResult;
class ContigMapData;
class OpticalMapData;
class AlignmentParams;

//Match the fragments from a single contig to the entire optical map
// return the best match as result
// if match_others is true, then put other match results in MatchResult
MatchResult * match(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     std::vector<MatchResult *> * pOthers, const AlignmentParams& alignParams);
MatchResult * match2(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     std::vector<MatchResult *> * pOthers, const AlignmentParams& alignParams);

// Light-weight alignment for permutation test.
// This alignment avoids filtering and avoids the construction of MatchResults.
double matchPermutationTest(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                        const AlignmentParams& alignParams);

//Match the fragments from a single contig to the entire optical map
// return the best match as result
// if match_others is true, then put other match results in MatchResult
MatchResult * matchLocal(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     std::vector<MatchResult *> * pOthers, const AlignmentParams& alignParams);

#endif
