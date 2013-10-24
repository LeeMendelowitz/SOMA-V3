#include "MapChunk.h"
#include "MapData.h"

using namespace std;

//size_t MapChunk::getStartIndex() const { return bFrag_ - map_->getFragsB(); }
//size_t MapChunk::getEndIndex() const { return eFrag_ - map_->getFragsB(); }

int sum(const FragDataConstIter bFrag, const FragDataConstIter eFrag)
{
    int s = 0;
    for(FragDataConstIter iter = bFrag; iter != eFrag; iter++)
        s += iter->size_;
    return s;
}

MapChunk::MapChunk(const MapData * map, const FragDataConstIter bFrag,
                                const FragDataConstIter eFrag) :
    map_(map),
    bFrag_(bFrag),
    eFrag_(eFrag),
    rank_(-1),
    size_(sum(bFrag, eFrag))
{}


MapChunk::MapChunk(const MapChunk& other) :
    map_(other.map_),
    bFrag_(other.bFrag_),
    eFrag_(other.eFrag_),
    next_(other.next_),
    prev_(other.prev_),
    rank_(other.rank_),
    size_(other.size_)
{}

bool MapChunk::isFirstQueryChunk() const { return getStartIndex() == 0; }
bool MapChunk::isLastQueryChunk() const { return getEndIndex() == map_->numFrags(); }
bool MapChunk::isBoundaryChunk() const { return isFirstQueryChunk() || isLastQueryChunk(); }

std::ostream& operator<<(std::ostream& o, const MapChunk& p)
{
    o << " Map: " << p.map_->getId()
      << " bInd: " << p.getStartIndex()
      << " eInd: " << p.getEndIndex()
      << " size (bp): " << p.size_;
    return o;
}
