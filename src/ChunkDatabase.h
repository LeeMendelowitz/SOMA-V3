/*
The ChunkDatabase stores reference fragments in a sorted array, which makes it possible
to efficiently find hits for a query.
*/
#ifndef FRAGDATABASE_H
#define FRAGDATABASE_H

#include <vector>
#include "OpticalMapData.h"
#include "seedTypes.h"

typedef std::pair<int, int> IntPair;
typedef std::vector<IntPair> IntPairVec;

class ChunkDatabase
{

    public:

    ChunkDatabase() {};
    ~ChunkDatabase();

    void addMap(const MapData * pMap, size_t maxInteriorMisses);
    void sortFrags();
    int  lowerBoundIndex(int q);
    bool getMapChunkHits(IntPairVec& query, std::vector<MapChunk*>& hits);
    bool getMapChunkHits(int lowerBound, int upperBound, std::vector<MapChunk*>& hits);

    MapChunkVec::const_iterator chunksB() const { return chunks_.begin();}
    MapChunkVec::const_iterator chunksE() const { return chunks_.end();}

    private:

    MapChunkVec chunks_;
};

#endif
