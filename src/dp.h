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
// TODO: Consider making classes to wrap the functions that create the scoreMatrices
// and compute the scores. This would avoid the need to constantly pass AlignmentParams
// around, or to consistently reference global variables.
class AlignmentParams
{
    public:
    AlignmentParams(double C_r_contig_in, double C_r_optical_in, double C_sigma_in,
                    double sigma2_in, int delta_in,
                    double smallFrag_in, double smallFragSlope_in,
                    double H_in, double T_in)
    {
        C_r_contig = C_r_contig_in;
        C_r_optical = C_r_optical_in;
        chi2Max = C_sigma_in * C_sigma_in;
        sigma2 = sigma2_in;
        delta = delta_in;

        smallFrag = smallFrag_in;
        smallFragSlope = smallFragSlope_in;

        // Shape parameters for parabola for local scoring function
        // Consider creating a class for 
        T2 = T_in * T_in;
        A = H_in / (T2); 

    }

    double C_r_contig; // Cost of missing a restriction site in the contig insilico map
    double C_r_optical; // Cost of missing a restriction site in the optical map
    double chi2Max; // The maximum allowed length difference for aligned fragments, in variance units.
    double sigma2; // The parameter for sigma^2 for computing the variance
    int delta; // Maximum number of unaligned sites inside an alignment block

    // Params for penalizing gapped contig fragments
    // or missed contig sites in small fragments
    double smallFrag;
    double smallFragSlope;

    // Params for local alignment
    double A;
    double T2;

};


//Match the fragments from a single contig to the entire optical map
// return the best match as result
// if match_others is true, then put other match results in MatchResult
MatchResult * match(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     vector<MatchResult *> * pOthers, bool forward, const AlignmentParams& alignParams);

// Light-weight alignment for permutation test.
// This alignment avoids filtering and avoids the construction of MatchResults.
double matchPermutationTest(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                        bool forward, const AlignmentParams& alignParams);

//Match the fragments from a single contig to the entire optical map
// return the best match as result
// if match_others is true, then put other match results in MatchResult
MatchResult * matchLocal(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     vector<MatchResult *> * pOthers, bool forward, const AlignmentParams& alignParams);


#endif
