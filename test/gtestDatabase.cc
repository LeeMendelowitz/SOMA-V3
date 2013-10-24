// Test the ChunkDatabase
#include "gtest/gtest.h"

#include <iostream>
using std::cout;
#include <vector>
using std::vector;

#include "ContigMapData.h"
#include "OpticalMapData.h"
#include "ChunkDatabase.h"
#include "MapChunkUtils.h"
#include "ScoreMatrixSeeded.h"
#include "seededDp.h"
#include "AlignmentParams.h"
#include "globals.h"

// To use a test fixture, derive a class from testing::Test.
class DatabaseTest : public testing::Test {
 protected:  // You should make the members protected s.t. they can be
             // accessed from sub-classes.

  // virtual void SetUp() will be called before each test is run.  You
  // should define it if you need to initialize the varaibles.
  // Otherwise, this can be skipped.
  virtual void SetUp()
  {
    readOpticalMaps("map100.opt", opticalMaps_);
    readContigMaps("contigs.silico", contigMaps_);
    std::cout << "Read " << contigMaps_.size() << " contig maps.\n";
    std::cout << "Read " << opticalMaps_.size() << " optical maps.\n";
    size_t maxInteriorMisses = 3;
    setMapChunks(opticalMaps_.front(), maxInteriorMisses);
    chunkDB_.addChunks(opticalMaps_.front()->getChunks());
    chunkDB_.sortFrags();
  }

  // virtual void TearDown() will be called after each test is run.
  // You should define it if there is cleanup work to do.  Otherwise,
  // you don't have to provide it.
  //
  virtual void TearDown() {
    for(size_t i = 0; i < opticalMaps_.size(); i++)
        delete opticalMaps_[i];
    opticalMaps_.clear();
    for(size_t i = 0; i < contigMaps_.size(); i++)
        delete contigMaps_[i];
    contigMaps_.clear();
  }

 
  void testFragConstructor()
  { }

  // Declares the variables your tests want to use.
  vector<ContigMapData*> contigMaps_;
  vector<OpticalMapData*> opticalMaps_;
  ChunkDatabase chunkDB_;
};


bool hitWithinBounds(const IntPairVec& bounds, MapChunk * hit)
{
    // NEEDS TO BE REWORKED

    /*
    FragDataVec::const_iterator d = hit->pFrag_;
    size_t i = 0;
    bool isMatch = true;
    for(; i < bounds.size(); i++, d++)
    {
        int lb = bounds[i].first;
        int ub = bounds[i].second;
        isMatch = isMatch && (d->size_ <= ub) && (d->size_ >= lb);
        if (!isMatch) break;
    }
    return isMatch;
    */
    return true;
}

bool hitWithinBounds(int lowerBound, int upperBound, MapChunk * hit)
{
    return (hit->size_ <= upperBound) && (hit->size_ >= lowerBound);
}

bool indexIsCorrect(MapChunkVec::const_iterator B, MapChunkVec::const_iterator E, int ind, int q)
{
    bool isOK = true;
    MapChunkVec::const_iterator iter = B + ind;

    // Check that this is the first element greater than q
    if (iter < E)
    {
        MapChunk * p = *iter;
        isOK = isOK &&  (p->size_ > q);
    }

    if (iter - 1 >= B)
    {
        MapChunk * p = *(iter-1);
        isOK = isOK && (p->size_ < q);
    }
    return isOK;
}


// When you have a test fixture, you define a test using TEST_F
// instead of TEST.

// Tests the default c'tor.
TEST_F(DatabaseTest, ConstructerSortTest) {
    MapChunkVec::const_iterator i = chunkDB_.chunksB();
    int lastSize = (*i)->size_;
    int lastRank = (*i)->rank_;
    ASSERT_EQ(0, lastRank);
    i++;
    for(; i != chunkDB_.chunksE(); i++)
    {
       MapChunk * pFrag = *i;
       ASSERT_GE(pFrag->size_, lastSize);
       ASSERT_EQ(lastRank + 1, pFrag->rank_);
       lastSize = pFrag->size_;
       lastRank++;
    }
}

template <class T>
bool testMembership(const vector<T>& vec, T test)
{
    for(typename vector<T>::const_iterator iter = vec.begin();
        iter != vec.end();
        iter++)
        if (*iter == test)
        {
            return true;
        }
    return false;
}

TEST_F(DatabaseTest, MapChunkPrevNextTest) {

    MapChunkVec::const_iterator i = chunkDB_.chunksB();
    for(; i != chunkDB_.chunksE(); i++)
    {
       MapChunk * thisChunk = *i;

       // Test that the next pointer is correct
       for(MapChunkVec::const_iterator iter = thisChunk->next_.begin();
           iter != thisChunk->next_.end();
           iter++)
       {
           MapChunk * nextChunk = *iter;
           ASSERT_EQ(nextChunk->bFrag_, thisChunk->eFrag_);
           ASSERT_TRUE(testMembership(nextChunk->prev_, thisChunk));
       }
    }
}

TEST_F(DatabaseTest, LowerBoundTest)
{
    std::vector<int> ints;
    ints.push_back(2000);
    ints.push_back(5000);
    ints.push_back(10000);
    ints.push_back(30000);
    ints.push_back(50000);
    ints.push_back(10000000);

    for (size_t i = 0; i < ints.size(); i++)
    {
        int q = ints[i];
        int ind = chunkDB_.lowerBoundIndex(q);
        bool isOK = indexIsCorrect(chunkDB_.chunksB(), chunkDB_.chunksE(), ind, q);
        ASSERT_TRUE(isOK);
    }

}

TEST_F(DatabaseTest, QueryTest1)
{
    // Make a query
    IntPairVec bounds;
    bounds.push_back(IntPair(30000,50000));
    bounds.push_back(IntPair(30000,50000));
    bounds.push_back(IntPair(30000,50000));

    vector<MapChunk*> hits;
    chunkDB_.getMapChunkHits(bounds, hits);

    cout << "Found " << hits.size() << " hits.\n";
    /*
    for(size_t i = 0; i < hits.size(); i++)
    {
        cout << *(hits[i]) << endl;
    }
    */

    for(size_t i=0; i < hits.size(); i++)
    {
        ASSERT_TRUE(hitWithinBounds(bounds, hits[i]));
    }

    cout << "Done with test." << endl;
}

TEST_F(DatabaseTest, QueryTest2)
{
    // Make a query
    vector<MapChunk*> hits;
    const int lb = 30000;
    const int ub = 50000;
    chunkDB_.getMapChunkHits(lb, ub, hits);

    cout << "Found " << hits.size() << " hits." << endl;
    /*
    for(size_t i = 0; i < hits.size(); i++)
    {
        cout << *(hits[i]) << "\n";
    }*/

    for(size_t i=0; i < hits.size(); i++)
    {
        ASSERT_TRUE(hitWithinBounds(lb, ub, hits[i]));
    }
}

TEST_F(DatabaseTest, CountCellsInPlay)
{
    float tol = 0.10;
    int minDelta = 1000;
    size_t maxInteriorMisses = 1;

    // Populate chunks for each query
    vector<MapChunkVec> mapChunks(contigMaps_.size());
    for (size_t i = 0; i < contigMaps_.size(); i++)
    {
        ContigMapData * cMap = contigMaps_[i];
        setMapChunks(cMap, maxInteriorMisses);
        RefToCoordSet refToCoordSet;
        calculateCellsInPlay(cMap->getChunks(), chunkDB_, tol, minDelta, refToCoordSet);
        ASSERT_TRUE(refToCoordSet.size() == 1);

        size_t numCells = 0;
        for(RefToCoordSet::const_iterator iter = refToCoordSet.begin(); iter != refToCoordSet.end(); iter++)
        {
            numCells += iter->second.size();
        }
        cout << "Got " << numCells << "cells.\n";
    }
}

TEST_F(DatabaseTest, GetScorePaths)
{
    float tol = 0.10;
    int minDelta = 1000;
    size_t maxInteriorMisses = 1;

    using seeded::ScoreMatrix;

    // Populate chunks for each query
    vector<MapChunkVec> mapChunks(contigMaps_.size());
    AlignmentParams alignParams(3, 5, 5, Constants::SIGMA, 3, 3, 0, 0, 0, 0);
    for (size_t i = 0; i < contigMaps_.size(); i++)
    {
        ContigMapData * cMap = contigMaps_[i];
        setMapChunks(cMap, maxInteriorMisses);

        RefToScorePathSteps refToScorePaths;  
        getScorePaths(cMap->getChunks(), chunkDB_, tol, minDelta, refToScorePaths);
        seeded::ScoreMatrix scoreMatrix(0,0);

        ASSERT_TRUE(refToScorePaths.size() == 1);

        size_t numScorePaths = 0;
        for(RefToScorePathSteps::iterator iter = refToScorePaths.begin(); iter != refToScorePaths.end(); iter++)
        {
            populateScoreMatrix(iter->second, cMap, iter->first, scoreMatrix);
            cout << "score matrix has size " << scoreMatrix.getSize() << ", capacity " << scoreMatrix.getCapacity() << endl;
            numScorePaths += iter->second.size();
            cout << "doing dp..." << flush;
            dp(scoreMatrix, alignParams);
            cout << "done.\n";
            cout << "number of filled cells: " << scoreMatrix.countFilledCells() << "\n";
            cout << "Max score: " << scoreMatrix.getMaxScore() << endl;
        }
        cout << "Got " << numScorePaths << "scorePaths.\n";
    }
}
