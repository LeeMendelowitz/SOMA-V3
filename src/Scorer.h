// Base Class to compute alignment scores for a chunk #ifndef SCORER_H
#ifndef SCORER_H
#define SCORER_H

#include <vector>
#include "mapTypes.h"
#include "AlignmentParams.h"

// Forward Declarations
class MatchResult;
class MatchedChunk;

class Score
{
    public:

    Score() :
        contig(0),
        optical(0),
        sizing(0)
    {}

    Score(double c, double o, double s) :
        contig(c),
        optical(o),
        sizing(s)
    {}
    double contig; // Score for contig (gap or site misses)
    double optical; // Score for optical (misses)
    double sizing; // Score for sizing error
};

class Scorer
{

    public:
    Scorer(const AlignmentParams& ap) : ap_(ap) {};
    virtual ~Scorer() {}

    virtual Score scoreGap(int gapSize) = 0;
    virtual Score scoreAlignment(const std::vector<FragData>::const_iterator cB,
                         const std::vector<FragData>::const_iterator cE,
                         const std::vector<FragData>::const_iterator oB,
                         const std::vector<FragData>::const_iterator oE,
                         bool boundaryFrag) = 0;

    Score scoreMatchedChunk(const MatchedChunk& chunk);

    // Score the matched chunks in a MatchResult
    void scoreMatchResult(MatchResult * pResult);


    protected:
        AlignmentParams ap_;
};

#endif
