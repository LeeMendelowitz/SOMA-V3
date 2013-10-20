#ifndef SEEDTYPES_H
#define SEEDTYPES_H

#include <vector>
#include <ostream>
#include "mapTypes.h"
#include "OpticalMapData.h"

// Class to represent a chunk from
// a query map or reference map.

class MapChunk;
typedef std::vector<FragData>::const_iterator FragDataConstIter;
typedef std::vector<MapChunk*> MapChunkVec;
typedef std::vector<MapChunk*>::iterator MapChunkVecIter;
typedef std::vector<MapChunk*>::const_iterator MapChunkVecConstIter;

class MapChunk
{
    public:

    // methods
    MapChunk(const MapData * map, const FragDataConstIter bFrag,
                                 const FragDataConstIter eFrag);
    MapChunk(const MapChunk& other);

    void setRank(int rank) { rank_ = rank; }
    int getRank() const { return rank_; }

    bool operator<(const MapChunk& other) const {
        return this->size_ < other.size_;
    }

    bool operator<(int fragSize) const {
        return this->size_ < fragSize;
    }

    size_t getStartIndex() const { return bFrag_ - map_->getFragsB(); }
    size_t getEndIndex() const { return eFrag_ - map_->getFragsB(); }

    // members
    const MapData * map_; // pointer to the map which owns the fragment
    const FragDataConstIter bFrag_; // pointer to the fragment
    const FragDataConstIter eFrag_; // pointer to the fragment
    MapChunkVec next_; // chunks which immediately follow this one in the map
    MapChunkVec prev_; // chunks which immediately preceed this on in the map
    int rank_; // rank of the fragment (sorted by size)
    const int size_; // bp in chunk
};


std::ostream& operator<<(std::ostream& o, const MapChunk& p);



#endif
