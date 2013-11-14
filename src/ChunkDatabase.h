/*
The ChunkDatabase stores reference fragments in a sorted array, which makes it possible
to efficiently find hits for a query.
*/
#ifndef CHUNKDATABASE_H
#define CHUNKDATABASE_H

#include <vector>
#include "OpticalMapData.h"
#include "MapChunk.h"

typedef std::pair<int, int> IntPair;
typedef std::vector<IntPair> IntPairVec;
typedef MapChunkVec::const_iterator MapChunkVecConstIter;
typedef std::pair<MapChunkVecConstIter, MapChunkVecConstIter> MapChunkVecConstIterPair;

class ChunkDatabase
{

    public:

    ChunkDatabase() {};
    ~ChunkDatabase();

    void addChunks(const MapChunkVec& chunks) { chunks_.insert(chunks_.end(), chunks.begin(), chunks.end()); }
    void sortFrags();
    int lowerBoundIndex(int q) const;
    MapChunkVecConstIter lowerBoundIter(int q) const;
    int upperBoundIndex(int q) const;
    int upperBoundIndex(int q, int lower) const;
    MapChunkVecConstIter upperBoundIter(int q) const;
    MapChunkVecConstIter upperBoundIter(int q, MapChunkVecConstIter lower) const;
    bool getMapChunkHits(IntPairVec& query, std::vector<MapChunk*>& hits) const;
    bool getMapChunkHits(int lowerBound, int upperBound, std::vector<MapChunk*>& hits) const;
    MapChunkVecConstIterPair getMapChunkHitCoords(int lowerBound, int upperBound) const;
    MapChunkVecConstIter chunksB() const { return chunks_.begin();}
    MapChunkVecConstIter chunksE() const { return chunks_.end();}

    private:

    MapChunkVec chunks_;
};

#endif
