// Abstract base class for MapData
#ifndef MAPDATA_H
#define MAPDATA_H

#include <vector>
#include <utility>
#include <iostream>
#include "mapTypes.h"
#include "MapChunk.h"
#include "MapChunkUtils.h"

typedef std::vector<FragData> FragDataVec;

class MapData
{
    public:
    MapData() {};
    MapData(const std::string& id) : id_(id) {};
    virtual ~MapData() { freeChunks(); };

    std::string getId() const { return id_; }

    // Get starting and ending bp for a fragment at the provided index.
    virtual int getStartBp(int ind, bool forward=true) const = 0;
    virtual int getEndBp(int ind, bool forward=true) const = 0;

    virtual int getLength() const = 0;

    // Return a constant reference to the map fragments
    virtual const FragDataVec& getFrags() const = 0;
    virtual FragDataVec::const_iterator getFragsB() const = 0;
    virtual FragDataVec::const_iterator getFragsE() const = 0;

    size_t numFrags() const { return getFrags().size(); }

    // Note: Chunk storage within MapData could be implemented better/safer
    void setChunks(MapChunkVec& chunks, MapChunkIndex& startToChunks, MapChunkIndex& endToChunks )
    {
        chunks_ = move(chunks);
        startToChunks_ = move(startToChunks);
        endToChunks_ = move(endToChunks);
    }

    void freeChunks() { for(auto c: chunks_) delete c; }
    const MapChunkVec& getChunksByStartIndex(int index) { return startToChunks_[index]; }
    const MapChunkVec& getChunksByEndIndex(int index) { return endToChunks_[index]; }
    const MapChunkVec& getChunks() { return chunks_; }

    private:
    void indexChunks();
    std::string id_;
    MapChunkVec chunks_;
    MapChunkIndex startToChunks_;
    MapChunkIndex endToChunks_;
};



#endif
