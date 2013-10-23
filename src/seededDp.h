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
#include "ScorePathStep.h"
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
}


typedef unordered_map<const MapData *, ScorePathStepVec> RefToScorePathSteps;
typedef RefToScorePathSteps::iterator RefToScorePathStepsIter;
typedef unordered_map<const MapData *, CoordSet> RefToCoordSet;
typedef RefToCoordSet::iterator RefToCoordSetIter;

void calculateCellsInPlay(const MapChunkVec& queryChunks, ChunkDatabase& chunkDB, float tol, int minDelta,
                          RefToCoordSet& refToCoordSet);
void getCells(const ScorePathStepVec& vec, CoordSet& coordSet);

#endif
