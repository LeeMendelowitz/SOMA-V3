#include <algorithm>
using std::sort;
using std::lower_bound;
using std::upper_bound;
#include <cassert>
#include <iostream>

#include "ChunkDatabase.h"

void incrementMapChunks(MapChunkVec& v);
void decrementMapChunks(MapChunkVec& v);
void filterNullMapChunks(MapChunkVec& v);
void filterMapChunksByRank(MapChunkVec&v, int lowerBound, int upperBound);
void filterMapChunksBySize(MapChunkVec&v, int lowerSize, int upperSize);

/*
Add fragments from the specified OpticalMapData to the 
database.
*/
void ChunkDatabase::addMap(const MapData * pMap, size_t maxInteriorMisses)
{

    // Create MapChunk's and put them in the chunks_ vector.

    typedef vector<MapChunkVec> OffsetVec;
    OffsetVec offset2Chunks(pMap->numFrags());

    MapChunkVec newChunks;
    newChunks.reserve(pMap->numFrags()*(maxInteriorMisses+1));

    const FragDataVec::const_iterator fragE = pMap->getFragsE();
    for (FragDataVec::const_iterator iter = pMap->getFragsB();
         iter != fragE;
         iter++)
    {
        for (size_t i = 0; i <= maxInteriorMisses; i++)
        {
            if (iter + i + 1 > fragE) break;
            MapChunk * f = new MapChunk(pMap, iter, iter + i + 1);
            assert(f->getStartIndex() < offset2Chunks.size());
            offset2Chunks[f->getStartIndex()].push_back(f);
            newChunks.push_back(f);
        }
    }

    // Set the neighbors 
    for(MapChunkVec::const_iterator citer = newChunks.begin();
        citer != newChunks.end();
        citer++)
    {
        MapChunk * pChunk = *citer;
        size_t offset = pChunk->getEndIndex();
        if (offset < offset2Chunks.size())
        {
            const MapChunkVec& successors = offset2Chunks[offset];
            for(MapChunkVec::const_iterator siter = successors.begin();
                siter != successors.end();
                siter++)
            {
                MapChunk * nextChunk = *siter;
                pChunk->next_.push_back(nextChunk);
                nextChunk->prev_.push_back(pChunk);
            }
        }
        chunks_.push_back(pChunk);
    }
}

ChunkDatabase::~ChunkDatabase()
{

    // Delete allocated memory
    for (MapChunkVec::iterator iter = chunks_.begin();
         iter != chunks_.end();
         iter++)
    {
        delete *iter;
    }

    chunks_.clear();

}

inline bool MapChunkCmp(const MapChunk * p1, const MapChunk* p2)
{
    return (*p1 < *p2);
}

inline bool MapChunkIntCmp(const MapChunk * p, int q)
{
    return (*p < q);
}

///////////////////////////////////////////////////////////////////////////////
/*
    sortFrags

    Sort the fragments, and assign the ranks.

*/
void ChunkDatabase::sortFrags()
{
    // Sort the fragments
    sort(chunks_.begin(), chunks_.end(), MapChunkCmp);

    // Assign ranks
    int rank = 0;
    for (MapChunkVec::iterator iter = chunks_.begin();
         iter != chunks_.end();
         iter++)
    {
        (*iter)->rank_ = rank++;
    }
}

inline int ChunkDatabase::lowerBoundIndex(int q)
{
    MapChunkVec::const_iterator iter = lower_bound(chunks_.begin(), chunks_.end(), q, MapChunkIntCmp);
    int index = iter - chunks_.begin();
    return index;
}

// Increment all MapChunks in the vector
void incrementMapChunks(MapChunkVec& v)
{
    const MapChunkVec::iterator vB = v.begin();
    const MapChunkVec::iterator vE = v.end();
    MapChunkVec newV;
    for(MapChunkVec::iterator i = vB; i != vE; i++)
    {
        MapChunk * c = *i;
        newV.insert(newV.end(), c->next_.begin(), c->next_.end());
    }
    filterNullMapChunks(newV);
    v.swap(newV);
}

// Decrement all MapChunks in the vector
void decrementMapChunks(MapChunkVec& v)
{
    // This needs to be fixed now that we allow mulitple successors
    /*
    for(MapChunkVec::iterator i = v.begin();
        i != v.end();
        i++)
    {
        MapChunk * old = *i;
        *i = old->pPrev_;
    }
    filterNullMapChunks(v);
    */
}

// Remove NULL MapChunks in the vector
void filterNullMapChunks(MapChunkVec& v)
{
    const MapChunk* pNull = NULL;
    MapChunkVec::iterator newEnd = remove(v.begin(), v.end(), pNull);
    size_t newSize = newEnd - v.begin();
    v.resize(newSize);
}

// Filter out MapChunk's that do not fall within the rank bounds.
void filterMapChunksByRank(MapChunkVec&v, int lowerBound, int upperBound)
{
    assert(lowerBound <= upperBound);
    MapChunkVec::iterator iter = v.begin();
    const MapChunkVec::iterator end = v.end();
    MapChunkVec::iterator loc = iter;
    while(iter != end)
    {
        MapChunk* p = *iter;
        if (p->rank_ >= lowerBound && p->rank_ < upperBound)
            *(loc++) = p;
        iter++;
    }

    size_t newSize = loc-v.begin();

    v.resize(newSize);
}

// Filter out MapChunk's that do not fall within the size bounds
void filterMapChunksBySize(MapChunkVec&v, int lowerSize, int upperSize)
{
    assert(lowerSize <= upperSize);
    MapChunkVec::iterator iter = v.begin();
    const MapChunkVec::iterator end = v.end();
    MapChunkVec::iterator loc = v.begin();
    while(iter != end)
    {
        MapChunk* p = *iter;
        if (p->size_ >= lowerSize && p->size_ < upperSize)
        {
            *loc = p;
            loc++;
        }
        iter++;
    }

    size_t newSize = loc - v.begin();
    v.resize(newSize);
}

bool ChunkDatabase::getMapChunkHits(IntPairVec& query, MapChunkVec& hits)
{

    if (query.empty())
        return false;

    hits.clear();
    const size_t numFrags = query.size();

    // Populate the hits vector using the first set of bounds.
    IntPairVec::const_iterator qi = query.begin();
    int lb = qi->first; // lower bound of query
    int ub = qi->second; // upper bound of query
    assert(lb <= ub);
    int lbi = lowerBoundIndex(lb);
    int ubi = lowerBoundIndex(ub);
    assert(lbi <= ubi);
    hits.reserve(ubi - lbi + 1);
    for(int i = lbi; i < ubi; i++)
        hits.push_back(chunks_[i]);

    qi++;
    for (; qi != query.end(); qi++)
    {

        // Increment all hits forward
        incrementMapChunks(hits);
        int lb = qi->first; // lower bound of query
        int ub = qi->second; // upper bound of query
        assert(lb <= ub);

        // Filter out MapChunks based on rank.
        filterMapChunksBySize(hits, lb, ub);

        if (hits.empty())
            break;
    }

    // Decrement the MapChunk's by number of fragments in the query.
    const size_t sizeBefore = hits.size();
    for (size_t i = 0; i < numFrags-1; i++)
        decrementMapChunks(hits);
    const size_t sizeAfter = hits.size();
    assert(sizeBefore == sizeAfter);

    return !hits.empty();
}

bool ChunkDatabase::getMapChunkHits(int lowerBound, int upperBound, MapChunkVec& hits)
{
    hits.clear();

    // Populate the hits vector using the first set of bounds.
    assert(lowerBound <= upperBound);
    int lbi = lowerBoundIndex(lowerBound);
    int ubi = lowerBoundIndex(upperBound);
    assert(lbi <= ubi);
    hits.reserve(ubi - lbi + 1);
    for(int i = lbi; i < ubi; i++)
        hits.push_back(chunks_[i]);
    return !hits.empty();
}
