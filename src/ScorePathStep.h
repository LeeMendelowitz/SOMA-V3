#ifndef SCOREPATHSTEP_H
#define SCOREPATHSTEP_H

#include "MapChunk.h"

class ScoreCell;

// ScorePathStep identifies one step through the dynamic programming table
// as the alignment of a queryChunk with a referenceChunk
class ScorePathStep
{
    public:

    ScorePathStep(const MapChunk * queryChunk, const MapChunk * refChunk) : 
        queryChunk_(queryChunk),
        refChunk_(refChunk),
        src_(nullptr),
        tgt_(nullptr) { };

    IntPair startCoord() const { return IntPair(queryChunk_->getStartIndex(), refChunk_->getStartIndex()); }
    IntPair endCoord() const { return IntPair(queryChunk_->getEndIndex(), refChunk_->getEndIndex()); }
    int queryStart() const { return queryChunk_->getStartIndex(); }
    int queryEnd() const { return queryChunk_->getEndIndex(); }
    int refStart() const { return refChunk_->getStartIndex(); }
    int refEnd() const { return refChunk_->getEndIndex(); }
    const MapData * getRefMap() const { return refChunk_->map_; }
    const MapData * getQueryMap() const { return queryChunk_->map_; }
    int getQuerySize() const { return queryChunk_->size_; }
    int getRefSize() const { return refChunk_->size_; }
    bool isFirstQueryChunk() const { return queryChunk_->getStartIndex() == 0; }
    bool isLastQueryChunk() const { return queryChunk_->getEndIndex() == queryChunk_->map_->numFrags(); }
    bool isBoundaryChunk() const { return isFirstQueryChunk() || isLastQueryChunk(); }
    void setSource(ScoreCell * src) { src_ = src;}
    void setTarget(ScoreCell * tgt) { tgt_ = tgt;}
    const MapChunk * getQueryChunk() { return queryChunk_; }
    const MapChunk * getRefChunk() { return refChunk_; }

    private:
    // This is a step from src_ to target_, where src_'s coordinates
    // are at the endCoords and tgt_ is at the startCoords. In other words,
    // this ScorePathStep represents a backpointer.
    const MapChunk * queryChunk_;
    const MapChunk * refChunk_;
    ScoreCell * src_;
    ScoreCell * tgt_;
};

typedef std::vector<ScorePathStep> ScorePathStepVec;

#endif
