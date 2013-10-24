#include <iostream>
#include <utility>
using std::move;
#include <limits>

#include "seededDp.h"
#include "ScoreMatrixSeeded.h"
#include "ScoreCell.h"
#include "scoringFunctions.h"

using namespace std;

const int MAX_INT = std::numeric_limits<int>::max();

//using seeded::ScoreMatrix;

// Return the cells in the dynamic programming table which are in play
// Place them in the refToCoordSet.
void calculateCellsInPlay(const MapChunkVec& queryChunks, ChunkDatabase& chunkDB, float tol, int minDelta,
                          RefToCoordSet& refToCoordSet)
{

    // Calculate lower and upper bounds for the queryChunks
    size_t numChunks = queryChunks.size();

    vector<int> lowerBounds(numChunks), upperBounds(numChunks);
    vector<MapChunkVecConstIterPair> matches(numChunks);
    for(size_t i = 0; i < numChunks; i++)
    {
        int chunkSize = queryChunks[i]->size_;
        int delta = (int) (tol*chunkSize + 0.5);
        if (delta < minDelta) delta = minDelta;
        lowerBounds[i] = chunkSize - delta;
        upperBounds[i] = chunkSize + delta;
        matches[i] = chunkDB.getMapChunkHitCoords(lowerBounds[i], upperBounds[i]);
    }

    // Collect the subpaths which are in play by reference map. 
    RefToScorePathSteps refToScorePathSteps;
    for(size_t i = 0; i < matches.size(); i++) {
        MapChunk * const queryChunk = queryChunks[i];
        const MapChunkVecConstIterPair& iterPair = matches[i];
        const MapChunkVecConstIter first = iterPair.first;
        const MapChunkVecConstIter last = iterPair.second;
        for(MapChunkVecConstIter iter = first; iter != last; iter++) {
            const MapChunk * refChunk = *iter;
            ScorePathStep sp(queryChunk, refChunk);
            RefToScorePathStepsIter spIter  = refToScorePathSteps.find(sp.getRefMap());
            if (spIter != refToScorePathSteps.end()) {
                spIter->second.push_back(sp);
            } else {
                ScorePathStepVec spVec(1, sp);
                refToScorePathSteps.insert(RefToScorePathSteps::value_type(sp.getRefMap(), spVec));
            }
        }
    }

    refToCoordSet.clear();

    // Collect the cells that are in play
    for(RefToScorePathSteps::const_iterator iter = refToScorePathSteps.begin();
        iter != refToScorePathSteps.end(); iter++)
    {
        CoordSet coordSet;
        getCells(iter->second, coordSet);
        refToCoordSet.insert(RefToCoordSet::value_type(iter->first, coordSet));
    }
}

// Compute the score paths using the chunkDB.
void getScorePaths(const MapChunkVec& queryChunks, ChunkDatabase& chunkDB, float tol, int minDelta, RefToScorePathSteps& refToScorePathSteps)
{

    // Calculate lower and upper bounds for the queryChunks
    size_t numChunks = queryChunks.size();

    vector<int> lowerBounds(numChunks), upperBounds(numChunks);
    vector<MapChunkVecConstIterPair> matches(numChunks);
    for(size_t i = 0; i < numChunks; i++)
    {
        int chunkSize = queryChunks[i]->size_;
        int delta = (int) (tol*chunkSize + 0.5);
        if (delta < minDelta) delta = minDelta;
        lowerBounds[i] = chunkSize - delta;
        upperBounds[i] = chunkSize + delta;
        // Boundary chunks are treat as arbitrarily truncated restriction fragments.
        // Set there upper limit to the maximum possible value.
        if (queryChunks[i]->isBoundaryChunk()) {
            upperBounds[i] = MAX_INT;
        }
        matches[i] = move(chunkDB.getMapChunkHitCoords(lowerBounds[i], upperBounds[i]));
    }

    // Collect the subpaths which are in play by reference map. 
    refToScorePathSteps.clear();
    for(size_t i = 0; i < matches.size(); i++) {
        MapChunk * const queryChunk = queryChunks[i];
        const MapChunkVecConstIterPair& iterPair = matches[i];
        const MapChunkVecConstIter first = iterPair.first;
        const MapChunkVecConstIter last = iterPair.second;
        for(MapChunkVecConstIter iter = first; iter != last; iter++) {
            const MapChunk * refChunk = *iter;
            ScorePathStep sp(queryChunk, refChunk);
            RefToScorePathStepsIter spIter  = refToScorePathSteps.find(refChunk->map_);
            if (spIter != refToScorePathSteps.end()) {
                spIter->second.emplace_back(queryChunk, refChunk);
            } else {
                ScorePathStep sp(queryChunk, refChunk);
                ScorePathStepVec spVec(1, move(sp));
                refToScorePathSteps.insert(RefToScorePathSteps::value_type(refChunk->map_, move(spVec)));
            }
        }
    }
}


// Populate the ScoreMatrix with the scorePathSteps.
// This is equivalent to setting the allowable edges in the DAG dynamic programming structure
void populateScoreMatrix(ScorePathStepVec& scorePathSteps, const MapData * queryMap, const MapData * refMap, seeded::ScoreMatrix& scoreMatrix)
{
    const size_t numQueryFrags = queryMap->numFrags();
    const size_t numRefFrags = refMap->numFrags();
    scoreMatrix.resize(numQueryFrags + 1, numRefFrags + 1);
    for(auto& scorePathStep : scorePathSteps)
    {
        scoreMatrix.addScorePathStep(scorePathStep);
    }
}

// Fill in scores/backpointers in the dynamic programming table.
void dp(seeded::ScoreMatrix& scoreMatrix, size_t maxMissesQuery, size_t maxMissesReference, const AlignmentParams& alignParams)
{
    const size_t numRows = scoreMatrix.getNumRows();
    const size_t numCols = scoreMatrix.getNumCols();

    // Initialize first row
    ScorePathStep defaultBackpointer(nullptr, nullptr);
    for (size_t j = 0; j < numCols; j++)
    {
        ScoreCell * pCell = scoreMatrix.getCell(j);
        pCell->score_ = 0.0;
        pCell->backPointer_ = defaultBackpointer;
    }

    // Initialize first column
    for (size_t row = 0; row < numRows; row++)
    {
        ScoreCell * pCell = scoreMatrix.getCell(row,0);
        pCell->score_ = -Constants::INF;
        pCell->backPointer_ = defaultBackpointer;
    }


    for (size_t i = 1; i < numRows; i++)
    {
        size_t rowOffset = i*numCols;
        size_t curIndex = rowOffset + 1;
        for(size_t j = 1; j < numCols; j++, curIndex++)
        {
            ScoreCell * pCell = scoreMatrix.getCell(curIndex);
            
            // Try all allowed backpointers.
            for (ScorePathStepVec::iterator stepIter = pCell->backPointers_.begin();
                 stepIter != pCell->backPointers_.end();
                 stepIter++)
            {
                const ScorePathStep& pathStep = *stepIter;
                ScoreCell * pTarget = pathStep.getTarget();
                if (pTarget->score_ == -Constants::INF)
                    continue;
                double chunkScore = scoringFunction(pathStep.getNumQueryFrags()-1, pathStep.getNumRefFrags()-1,
                                             pathStep.getQuerySize(), pathStep.getRefSize(),
                                             pathStep.isBoundaryChunk(), alignParams);
                double newScore = pTarget->score_ + chunkScore;
                if (newScore > pCell->score_)
                {
                    pCell->score_ = newScore;
                    pCell->backPointer_ = pathStep;
                }
            }
        }
    }
}

/*
class ScoreTableTracker
{
    public:

    ScoreTableTracker(const MapData * queryMap, const MapData * refMap, 
                      size_t maxInteriorMissedSites) :
        queryMap_(queryMap),
        refMap_(refMap),
        maxInteriorMissedSites_(maxInteriorMissedSites)
    { };

    void processScorePathSteps(ScorePathStepVec& scorePaths);

    // Get or create a score cell.
    ScoreCell * getScoreCell(const IntPair& coord)
    {
        ScoreCellMap::iterator iter = scoreCellMap_.find(coord);
        if(iter == scoreCellMap_.end())
        {
            ScoreCell * cell = new ScoreCell(coord);
            scoreCells_.push_back(cell);
            scoreCellMap_.insert(ScoreCellMap::value_type(coord, cell));
            return cell;
        }
        return iter->second;
    }

    private:
    const MapData * queryMap_;
    const MapData * refMap_;
    ScoreCellMap scoreCellMap_;
    ScoreCellVec scoreCells_;
    size_t maxInteriorMissedSites_;
};

// Process score paths for a single query & single reference
// to identify cells that are in play.
void ScoreTableTracker::processScorePathSteps(ScorePathStepVec& scorePaths)
{
    for (auto sp : scorePaths)
    {
        ScoreCell * startCell = getScoreCell(sp.startCoord());
        ScoreCell * endCell = getScoreCell(sp.endCoord());
        startCell->forwardPointers_.push_back(endCell);
        endCell->forwardPointers_.push_back(startCell);
    }
    
}
*/

// Extract cells that are in play for a vector of ScorePathSteps.
void getCells(const ScorePathStepVec& vec, CoordSet& coordSet)
{
    for(size_t i = 0; i < vec.size(); i++)
    {
        const ScorePathStep& sp = vec[i];
        coordSet.insert(sp.endCoord());
    }
}
