#ifndef LOCAL_ALIGNMENT
#define LOCAL_ALIGNMENT

// Scoring functions for local alignment
// nContigSites: number of unaligned contig sites
// nOpticalSites: number of unaligned optical sites
double localScoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag)
{
    // Model assumes that Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size 
    double delta = (opticalLength - contigLength);
    double var = contigLength*Constants::SIGMA2;
    double C = 10.0;
    double A = 0.5;
    double sizeScore = A*((delta*delta/var) - C);
    double missPenalty = nContigSites*2.0 + nOpticalSites*2.0;
    return sizeScore + missPenalty;
}


#endif
