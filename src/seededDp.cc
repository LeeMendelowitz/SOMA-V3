#include <iostream>

#include "seededDp.h"
#include "ScoreCell.h"

using namespace std;

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
