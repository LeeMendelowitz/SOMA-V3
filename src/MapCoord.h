// Defines an interval in an Map

#ifndef MAPCOORD_H
#define MAPCOORD_H

#include "MapData.h"
#include "exception.h"

#include <sstream>

class MapCoord
{
    public:

    MapCoord() :
        pMap_(NULL),
        start_(-1),
        end_(-1)
    {}

    MapCoord(const MapData * pMap, int start, int end) :
        pMap_(pMap),
        start_(start),
        end_(end)
    {}

    int getLengthBp() const
    {
        return pMap_->getEndBp(end_-1) - pMap_->getStartBp(start_);
    }

    int getNumFrags() const
    {
        return end_ - start_;
    }

    int getStart() const
    {
        return start_;
    }

    int getEnd() const
    {
        return end_;
    }

    std::vector<FragData>::const_iterator getB() const
    {
        const std::vector<FragData>& frags = pMap_->getFrags();
        const size_t N = frags.size();
        if ((start_ < 0) || ((size_t) start_ > N))
        {
            std::ostringstream oss;
            oss << "Invalid MapCoord: "
                << "map: " << pMap_->getId() << " "
                << "numFrags: " << N
                << " start_: " << start_
                << " end_: " << end_;
            throw( Exception(oss.str()) );
        }

        return frags.begin() + start_;
        //return pMap_->getFrags().begin() + start_;
    }

    std::vector<FragData>::const_iterator getE() const
    {
        const std::vector<FragData>& frags = pMap_->getFrags();
        const size_t N = frags.size();
        if ((end_ < 0) || ((size_t) end_ > N))
        {
            std::ostringstream oss;
            oss << "Invalid MapCoord: "
                << "map: " << pMap_->getId() << " "
                << "numFrags: " << N
                << " start_: " << start_
                << " end_: " << end_;
            throw( Exception(oss.str()) );
        }
        return frags.begin() + end_;
        //return pMap_->getFrags().begin() + end_;
    }

    const MapData * getMap() const { return pMap_; }

    const MapData * pMap_;
    int start_; // start index, inclusive
    int end_; // end index, exclusive
};

#endif
