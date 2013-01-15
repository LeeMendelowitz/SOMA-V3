// Read input maps from a map file
#ifndef MAPREADER_H
#define MAPREADER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "mapTypes.h"

struct MapInput
{
    std::vector<FragData> frags_;
    std::string mapId_;
    int length_;
    int numFrags_;

    void reset();
};

class MapReader
{
    public:

    MapReader(const std::string& fileName);
    ~MapReader();

    bool next(MapInput& data);
    bool read(std::vector<MapInput>& data);

    private:
    void readLine(const std::string& line, MapInput& data);

    std::string fileName_;
    std::ifstream ifs_;
};

#endif
