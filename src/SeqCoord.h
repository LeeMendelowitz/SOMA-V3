// Defines an interval of a map, defined by bp positions
#ifndef SEQCOORD_H
#define SEQCOORD_H

#include "MapData.h"

class SeqCoord
{
    public:
    SeqCoord() :
        pMap_(NULL),
        start_(-1),
        end_(-1)
    {}

    SeqCoord(const MapData * pMap, int start, int end) :
        pMap_(pMap),
        start_(start),
        end_(end)
    {}

    int getLength() const 
    {
        return end_ - start_;
    }

    int getStartBp() const
    {
        return start_;
    }

    int getEndBp() const
    {
        return end_;
    }

    const MapData * getMap() const
    {
        return pMap_;
    }

    const MapData * pMap_;
    int start_; // inclusive
    int end_; // exclusive

};

#endif
