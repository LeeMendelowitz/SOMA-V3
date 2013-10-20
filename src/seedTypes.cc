#include "seedTypes.h"

using namespace std;

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
{};


MapChunk::MapChunk(const MapChunk& other) :
    map_(other.map_),
    bFrag_(other.bFrag_),
    eFrag_(other.eFrag_),
    next_(other.next_),
    prev_(other.prev_),
    rank_(other.rank_),
    size_(other.size_)
{};

std::ostream& operator<<(std::ostream& o, const MapChunk& p)
{
    o << " Map: " << p.map_->getId()
      << " bInd: " << p.getStartIndex()
      << " eInd: " << p.getEndIndex()
      << " size (bp): " << p.size_;
    return o;
}
