// dp.h
#ifndef DP_H
#define DP_H

#include <sstream>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <sstream>
#include <utility>

#include "exception.h"
#include "globals.h"
#include "utils.h"
#include "match.h"
#include "array2d.h"
#include "MatchResult.h"


// Parameters used for an alignment
// This is an input to the match function
class AlignmentParams
{
    public:
    double C_s; // Cost of editing a single base
    double C_r_contig; // Cost of missing a restriction site in the contig insilico map
    double C_r_optical; // Cost of missing a restriction site in the optical map
    double C_sigma; // The maximum allowed length difference for aligned fragments, in standard deviation units
    int delta; // Maximum number of unaligned sites inside an alignment block

    AlignmentParams(double C_s_in, double C_r_contig_in, double C_r_optical_in, double C_sigma_in, int delta_in)
    {
        C_s = C_s_in;
        C_r_contig = C_r_contig_in;
        C_r_optical = C_r_optical_in;
        C_sigma = C_sigma_in;
        delta = delta_in;
    }
};



//Match the fragments from a single contig to the entire optical map
// return the best match as result
// if match_others is true, then put other match results in MatchResult
MatchResult * match(const ContigMapData * pContigMapData, const OpticalMapData * pOpticalMapData,
                     vector<MatchResult *> * pOthers, bool forward, const AlignmentParams& alignParams);

// Light-weight alignment for permutation test.
// This alignment avoids filtering and avoids the construction of MatchResults.
double matchPermutationTest(const ContigMapData * pContigMapData, const OpticalMapData * pOpticalMapData,
                        bool forward, const AlignmentParams& alignParams);


#endif
