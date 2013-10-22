#ifndef SEEDEDDP_H
#define SEEDEDDP_H

#include <vector>
using std::vector;
#include <unordered_map>
using std::unordered_map;
#include <unordered_set>
using std::unordered_set;
#include <functional>

#include "mapTypes.h"
#include "MapChunk.h"
#include "ChunkDatabase.h"

typedef std::pair<int, int> IntPair;
typedef std::unordered_set<IntPair> CoordSet;

// Hash function declarations
namespace std {

    template <>
    struct hash< IntPair > {
        public:
        size_t operator()(const IntPair& x) const throw() {
             std::hash<int> hasher;
             return (hasher(x.first << 16) ^ hasher(x.second));
        }
    };
    /*
    template <>
    struct hash< ScoreCell > {
        public:
        size_t operator()(const ScoreCell& x) const throw();
    };
    */
}

// ScorePath identifies one step through the dynamic programming table
// as the alignment of a queryChunk with a referenceChunk
class ScorePath
{
    public:

    ScorePath(const MapChunk * queryChunk, const MapChunk * refChunk) : 
        queryChunk_(queryChunk),
        refChunk_(refChunk) { };

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

    const MapChunk * queryChunk_;
    const MapChunk * refChunk_;
};

// ScoreCell identifies one cell in the dynamic programming table and
// identifies its successors/predecessors.
class ScoreCell;
typedef std::vector<ScoreCell *> ScoreCellVec;
class ScoreCell
{
    public:
    ScoreCell(int q, int r) :
        q_(q), r_(r) { };
    ScoreCell(IntPair ip) :
        q_(ip.first), r_(ip.second) { };

    IntPair key() const { return IntPair(q_, r_);}
    bool operator==(const ScoreCell& other) const { return (key() == other.key()); } 
    size_t inDegree() const { return backPointers_.size(); };
    size_t outDegree() const { return forwardPointers_.size(); };

    const int q_;
    const int r_;
    ScoreCellVec backPointers_;
    ScoreCellVec forwardPointers_;
};

/*
template <>
size_t std::hash<ScoreCell>::operator()(const ScoreCell& x) const throw() {
     std::hash<IntPair> hasher;
     return hasher(x.key());
}
*/

typedef unordered_map<IntPair, ScoreCell *> ScoreCellMap;

typedef std::vector<ScorePath> ScorePathVec;
typedef unordered_map<const MapData *, ScorePathVec> RefToScorePaths;
typedef RefToScorePaths::iterator RefToScorePathsIter;
typedef unordered_map<const MapData *, CoordSet> RefToCoordSet;
typedef RefToCoordSet::iterator RefToCoordSetIter;

void calculateCellsInPlay(const MapChunkVec& queryChunks, ChunkDatabase& chunkDB, float tol, int minDelta,
                          RefToCoordSet& refToCoordSet);
void getCells(const ScorePathVec& vec, CoordSet& coordSet);

#endif
