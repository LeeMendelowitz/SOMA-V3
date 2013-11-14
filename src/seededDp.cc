#include <iostream>
#include <utility>
using std::move;
#include <limits>
#include <algorithm>
using std::min;

#include "MatchResult.h"
#include "seededDp.h"
#include "ScoreMatrixSeeded.h"
#include "StandardMatchMaker.h"
#include "ScoreCell.h"
#include "scoringFunctions.h"
#include "GlobalScorer.h"
#include "Clock.h"

#define SEEDED_DP_DEBUG 1

using seeded::ScoreMatrix;

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
void getScorePaths(const MapChunkVec& queryChunks, const ChunkDatabase& chunkDB, float tol, int minDelta, RefToScorePathSteps& refToScorePathSteps)
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
void setScoreMatrixPathSteps(ScorePathStepVec& scorePathSteps, const MapData * queryMap, const MapData * refMap, seeded::ScoreMatrix& scoreMatrix)
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
void fillScoreMatrix(seeded::ScoreMatrix& scoreMatrix, const AlignmentParams& alignParams)
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
    for (size_t row = 1; row < numRows; row++)
    {
        ScoreCell * pCell = scoreMatrix.getCell(row,0);
        pCell->score_ = -Constants::INF;
        pCell->backPointer_ = defaultBackpointer;
    }

    // Initialize rest
    for (size_t row = 1; row < numRows; row++) {
        for(size_t col = 1; col < numCols; col++) {
            ScoreCell * pCell = scoreMatrix.getCell(row,col);
            pCell->score_ = -Constants::INF;
            pCell->backPointer_ = defaultBackpointer;
        }
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
                assert(pathStep.getSource() == pCell);
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

///////////////////////////////////////////////////////////////////////////////////////////////
// Match the fragments from a single contig to all of the reference maps.
// return the best match as a new MatchResult
// if pAllMatches is not NULL, add all matches (inlcuding the best one) to pAllMatches vector.
MatchResult *  seededMatch(ContigMapData * cMap, const ChunkDatabase& chunkDB, 
                     ScoreMatrix& scoreMatrix, vector<MatchResult *> * pAllMatches, const AlignmentParams& alignParams)
{

    Clock clock;
    MatchResult * bestMatch = nullptr;
    bool contigIsForward = cMap->isForward();
    const vector<FragData>& contigFrags = cMap->getFrags();

    // Extract chunks from the query
    setMapChunks(cMap, alignParams.maxChunkMissesQuery);

    // Get compatible alignments between query and reference

    RefToScorePathSteps refToScorePaths;
    clock.lap();
    getScorePaths(cMap->getChunks(), chunkDB, alignParams.tol, alignParams.minDelta, refToScorePaths);
    std::cout << "Getting Score Paths: " << clock.lap() << " elapsed.\n";

    GlobalScorer scorer(alignParams);
    StandardMatchMaker matchMaker(opt::maxMatchesPerContig, opt::minContigHits,
                                  opt::minLengthRatio, opt::maxMissRateContig, opt::avgChi2Threshold);
    size_t numScorePaths = 0;

    // Build matches from the score matrix
    MatchResultPtrVec allMatches;

    clock.lap();
    for(RefToScorePathSteps::iterator iter = refToScorePaths.begin(); iter != refToScorePaths.end(); iter++)
    {
        auto refMap = iter->first;
        auto& scorePathSteps = iter->second;

        Clock clock;

        // Reset the score matrix
        clock.lap();
        scoreMatrix.resetCells();
        std::cout << "Reset: " << clock.lap() << "\n";
    

        ////////////////////////////////////////
        // Perform dynamic programming
        clock.lap();
        setScoreMatrixPathSteps(scorePathSteps, cMap, refMap, scoreMatrix); // set edges (forward pointers / back pointers)
        std::cout << "Set: " << clock.lap() << "\n";

        clock.lap();
        fillScoreMatrix(scoreMatrix, alignParams); // traverse table, and fill in the best score and best back pointer
        std::cout << "Fill: " << clock.lap() << "\n";

        clock.lap();
        MatchResultPtrVec matches;
        matchMaker.makeMatches(scoreMatrix, scorer, matches, refMap, cMap, contigIsForward);
        allMatches.insert(allMatches.end(), matches.begin(), matches.end());
        std::cout << "Made " << matches.size() << " matches."
                  << " Elapsed: " << clock.lap() << "\n";
    }

    
    std::cout << "-------------------------------------\n"
              << "Finished making matches. Made  "
              << allMatches.size() << " matches. "
              << clock.lap() << " elapsed.\n";

    clock.lap();
    sort(allMatches.begin(), allMatches.end(), MatchResult::compareScore);
    reverse(allMatches.begin(), allMatches.end());

    // Set the matches to be returned
    if (!allMatches.empty())
    {
        // The first match is the best match
        bestMatch = allMatches.front();

        MatchResultPtrVec::iterator E = allMatches.begin()+1;
        if (pAllMatches)
        {
            int offset = min((int) allMatches.size(), opt::maxMatchesPerContig);
            E = allMatches.begin() + offset;
            pAllMatches->insert(pAllMatches->end(), allMatches.begin(), E);
        }

        // Delete all other matches other than those that are returned
        for(auto iter = E; iter != allMatches.end(); iter++)
        {
            delete *iter;
        }
    }
    std::cout << "Sort and delete matches: " << clock.lap() << " seconds.\n";

    return bestMatch;
}
