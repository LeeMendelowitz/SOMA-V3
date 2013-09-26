#ifndef SEEDTYPES_H
#define SEEDTYPES_H

#include <vector>

#include "mapTypes.h"
#include "OpticalMapData.h"

class FragPtr
{
    public:

    // methods
    FragPtr(const OpticalMapData * map, const vector<FragData>::const_iterator pFrag);
    FragPtr(const FragPtr& other);

    void setRank(int rank) { rank_ = rank; }
    int getRank() const { return rank_; }

    bool operator<(const FragPtr& other) const {
        return this->pFrag_->size_ < other.pFrag_->size_;
    }

    // members
    std::vector<FragData>::const_iterator pFrag_; // pointer to the fragment
    FragPtr * pNext_;
    FragPtr * pPrev_;
    const OpticalMapData * map_; // pointer to the map which owns the fragment
    int rank_; // rank of the fragment (sorted by size)
};






#endif
