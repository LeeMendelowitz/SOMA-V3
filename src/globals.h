#ifndef GLOBALS_H
#define GLOBALS_H
#include <string>
#include <map>
#include <vector>

using namespace std;

namespace Constants
{
    extern const double INF;
    extern const int MAX_GAP_SIZE;
    extern const double MAX_LOST_FRAC; // maximum fraction of contig that is allowed to be lost
    extern const int FRAG_CUTOFF; // small fragment cutoff
    extern const double SIGMA2; // Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size. Units: (bp)
    extern const double SIGMA; // sqrt of SIGMA2
}

// Options set by ParseArgs
namespace Options
{
    extern string silicoMap;
    extern vector<string> opticalMapList;
    extern bool   circular;
//    extern bool allowFalseCuts;
//    extern bool   matchFragmentOnce;
    extern bool noReverse;
    extern bool oneToOneMatch;
//    extern bool allowGaps;
    extern double pThreshold;
    extern string outputPrefix;
    extern double sdMin;
    extern double sdMax;
    extern double C_r_optical;
    extern double C_r_contig;
    //extern double C_s;
    extern double maxMissRateContig;
    extern double minLengthRatio;
    extern double avgChi2Threshold;
    extern int maxMatchesPerContig;
//    extern double falseCutRate;
//    extern int maxGapSize;
    extern int numPermutationTrials;
    extern int numThreads;
    extern int delta; // maximum number of unaligned sites in interior of alignment block
    extern int smallFrag; // Cutoff of small fragments (bp)
    extern double smallFragSlope;
}

namespace Globals
{
    typedef std::map<char, char> DNAMap; // map from a base to its complement

    // Global Variables set by main
    extern DNAMap * pDNAMap;
    
    // initialize global variables
    void initialize();

    // Free dynamically allocated global variables
    void free();
}

#endif
