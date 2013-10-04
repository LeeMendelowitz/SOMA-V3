#include "seedTypes.h"

FragPtr::FragPtr(const OpticalMapData * map, const vector<FragData>::const_iterator pFrag) :
    pFrag_(pFrag),
    map_(map),
    pNext_(0),
    pPrev_(0),
    rank_(-1)
{};


FragPtr::FragPtr(const FragPtr& other) :
    pFrag_(other.pFrag_),
    map_(other.map_),
    pNext_(other.pNext_),
    pPrev_(other.pPrev_),
    rank_(other.rank_)
{};

std::ostream& operator<<(std::ostream& o, const FragPtr& p)
{
    o << " Map: " << p.map_->getId()
      << " index: " << p.getIndex()
      << " size (bp): " << p.pFrag_->size_;
    return o;
}