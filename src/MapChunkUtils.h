#ifndef CHUNKUTILS_H
#define CHUNKUTILS_H

#include <cstdlib>

class MapData;
void setMapChunks(MapData * pMap, size_t maxInteriorMisses = 0);

/*
class MapChunkGenerator
{
    public:

    MapChunkGenerator(const MapData * pMap, size_t maxInteriorMisses = 0);

    ~MapChunkGenerator();

    MapChunkVecConstIter getChunksB() { return chunks_.begin(); }
    MapChunkVecConstIter getChunksE() { return chunks_.end(); }
    MapChunkVecBounds getChunks() {
        return MapChunkVecBounds(chunks_.begin(), chunks_.end());
    }
    MapChunkVecBounds getChunksByStart(int index) {
        // Check if index is within bounds?
        const MapChunkVec& chunkVec = startToChunks_[index];
        return MapChunkVecBounds(chunkVec.begin(), chunkVec.end());
    }
    MapChunkVecBounds getChunksByEnd(int index) {
        // Check if index is within bounds?
        const MapChunkVec& chunkVec = endToChunks_[index];
        return MapChunkVecBounds(chunkVec.begin(), chunkVec.end());
    }

    private:
    void genChunks();

    const MapData * map_;
    size_t maxInteriorMisses_;
    IndexToChunks startToChunks_;
    IndexToChunks endToChunks_;
    MapChunkVec chunks_;
};
*/

#endif
