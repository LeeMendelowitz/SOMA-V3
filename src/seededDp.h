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
namespace seeded {class ScoreMatrix;}
class AlignmentParams;
class MatchResult;

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
void getScorePaths(const MapChunkVec& queryChunks, const ChunkDatabase& chunkDB, float tol, int minDelta, RefToScorePathSteps& refToScorePathSteps);
void getCells(const ScorePathStepVec& vec, CoordSet& coordSet);
void setScoreMatrixPathSteps(ScorePathStepVec& scorePathSteps, const MapData * queryMap, const MapData * refMap, seeded::ScoreMatrix& scoreMatrix);
void fillScoreMatrix(seeded::ScoreMatrix& scoreMatrix, const AlignmentParams& alignParams);
MatchResult * seededMatch(ContigMapData * cMap, const ChunkDatabase& chunkDB, seeded::ScoreMatrix& scoreMatrix,
                          vector<MatchResult *> * pAllMatches, const AlignmentParams& alignParams);

#endif
