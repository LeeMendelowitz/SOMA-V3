// Abstract base class for MapData
#ifndef MAPDATA_H
#define MAPDATA_H

#include <vector>
#include "mapTypes.h"


typedef std::vector<FragData> FragDataVec;

class MapData
{
    public:
    MapData() {};
    MapData(const std::string& id) : id_(id) {};
    virtual ~MapData() {};

    std::string getId() const { return id_; }

    // Get starting and ending bp for a fragment at the provided index.
    virtual int getStartBp(int ind, bool forward=true) const = 0;
    virtual int getEndBp(int ind, bool forward=true) const = 0;

    virtual int getLength() const = 0;

    // Return a constant reference to the map fragments
    virtual const FragDataVec& getFrags() const = 0;

    private:
    std::string id_;
};


#endif
