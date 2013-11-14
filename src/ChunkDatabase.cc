#include <algorithm>
using std::sort;
using std::lower_bound;
using std::upper_bound;
#include <cassert>
#include <iostream>

#include "ChunkDatabase.h"
#include "MapChunkUtils.h"

void incrementMapChunks(MapChunkVec& v);
void decrementMapChunks(MapChunkVec& v);
void filterNullMapChunks(MapChunkVec& v);
void filterMapChunksByRank(MapChunkVec&v, int lowerBound, int upperBound);
void filterMapChunksBySize(MapChunkVec&v, int lowerSize, int upperSize);

/*
Add fragments from the specified OpticalMapData to the 
database.
*/

ChunkDatabase::~ChunkDatabase()
{

    // Delete allocated memory
    for (MapChunkVec::iterator iter = chunks_.begin();
         iter != chunks_.end();
         iter++)
    {
//        delete *iter;
    }

    chunks_.clear();

}

inline bool MapChunkCmp(const MapChunk * p1, const MapChunk* p2)
{
    return (*p1 < *p2);
}

inline bool MapChunkIntCmp(const MapChunk * p, int q)
{
    return p->size_ < q;
}
inline bool IntMapChunkCmp(int q, const MapChunk * p)
{
    return q < p->size_;
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

inline int ChunkDatabase::lowerBoundIndex(int q) const
{
    MapChunkVecConstIter iter = lower_bound(chunks_.begin(), chunks_.end(), q, MapChunkIntCmp);
    int index = iter - chunks_.begin();
    return index;
}

inline MapChunkVecConstIter ChunkDatabase::lowerBoundIter(int q) const
{
    return lower_bound(chunks_.begin(), chunks_.end(), q, MapChunkIntCmp);
}

inline int ChunkDatabase::upperBoundIndex(int q) const
{
    MapChunkVecConstIter iter = upper_bound(chunks_.begin(), chunks_.end(), q, IntMapChunkCmp);
    int index = iter - chunks_.begin();
    return index;
}

inline int ChunkDatabase::upperBoundIndex(int q, int lowerIndex) const
{
    MapChunkVecConstIter iter = upper_bound(chunks_.begin() + lowerIndex, chunks_.end(), q, IntMapChunkCmp);
    int index = iter - chunks_.begin();
    return index;
}

inline MapChunkVecConstIter ChunkDatabase::upperBoundIter(int q) const
{
    return upper_bound(chunks_.begin(), chunks_.end(), q, IntMapChunkCmp);
}

inline MapChunkVecConstIter ChunkDatabase::upperBoundIter(int q, MapChunkVecConstIter lower) const
{
    MapChunkVecConstIter chunksE = chunks_.end();
    //return upper_bound(chunks_.begin(), chunks_.end(), q, MapChunkIntCmp);
    return upper_bound(lower, chunksE, q, IntMapChunkCmp);
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

bool ChunkDatabase::getMapChunkHits(IntPairVec& query, MapChunkVec& hits) const
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
    int ubi = upperBoundIndex(ub, lbi);
    assert(lbi <= ubi);
    hits.reserve(ubi - lbi);
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

bool ChunkDatabase::getMapChunkHits(int lowerBound, int upperBound, MapChunkVec& hits) const
{
    hits.clear();

    // Populate the hits vector using the first set of bounds.
    assert(lowerBound <= upperBound);
    int lbi = lowerBoundIndex(lowerBound);
    int ubi = upperBoundIndex(upperBound, lbi);
    assert(lbi <= ubi);
    hits.reserve(ubi - lbi);
    for(int i = lbi; i < ubi; i++)
        hits.push_back(chunks_[i]);
    return !hits.empty();
}

MapChunkVecConstIterPair ChunkDatabase::getMapChunkHitCoords(int lowerBound, int upperBound) const
{
    MapChunkVecConstIter lower = lowerBoundIter(lowerBound);
    MapChunkVecConstIter upper = upperBoundIter(upperBound);
    return MapChunkVecConstIterPair(lower, upper);
}
