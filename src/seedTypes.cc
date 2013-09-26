#include "seedTypes.h"

FragPtr::FragPtr(const OpticalMapData * map, const vector<FragData>::const_iterator pFrag) :
    map_(map),
    pFrag_(pFrag),
    pNext_(0),
    pPrev_(0),
    rank_(-1)
{};


FragPtr::FragPtr(const FragPtr& other) :
    pFrag_(other.pFrag_),
    pNext_(other.pNext_),
    pPrev_(other.pPrev_),
    map_(other.map_),
    rank_(other.rank_)
{};



