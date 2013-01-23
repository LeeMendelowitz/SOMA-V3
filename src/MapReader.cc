#include "MapReader.h"
#include "utils.h"
#include "exception.h"

using namespace std;

#define MAPREADER_DEBUG 0

void MapInput::reset()
{
    length_ = 0;
    numFrags_ = 0;
    mapId_ = "";
    frags_.clear();
}

std::ostream& operator<<(std::ostream& os, const MapInput& map)
{
    os << map.mapId_ << " "
       << map.numFrags_ << " "
       << map.length_ << ". ";
    for(std::vector<FragData>::const_iterator i = map.frags_.begin();
        i != map.frags_.end();
        i++)
    {
        os << i->size_ << " ";
    }
    return os;
}



MapReader::MapReader(const std::string& fileName) : fileName_(fileName) {
    ifs_.open(fileName_.c_str(), ios_base::in);
    if (ifs_.fail())
    {
        ostringstream msg;
        msg << "Error: Could not open map file " << fileName_ << ".";
        throw(Exception(msg.str()));
        return;
    }
}

MapReader::~MapReader()
{
    ifs_.close();
}

bool MapReader::next(MapInput& data)
{
    string line;
    bool success = getline(ifs_, line);
    if (!success) return false;
    readLine(line, data);
    return true;
}

bool MapReader::read(std::vector<MapInput>& dataVec)
{
    dataVec.clear();

    MapInput data;
    while(next(data))
    {
        dataVec.push_back(data);
    }
    return !dataVec.empty();
}

// parse an input line from a map file
void MapReader::readLine(const std::string& line, MapInput& data)
{
    const char delim = '\t';
    const size_t minFields = 3;

    data.reset();

    vector<string> fields;
    splitString(line, delim, &fields);

    bool lineOK = true;

    if (fields.size() < minFields) lineOK = false;

    vector<string>::const_iterator iter = fields.begin();
    const vector<string>::const_iterator E = fields.end();

    istringstream iss;
    if (lineOK)
    {
        // Read map id
        data.mapId_ = *iter++;

        // Read map length
        iss.clear();
        iss.str(*iter++);
        iss >> data.length_;
        lineOK = lineOK && !iss.fail();

        // Read number of frags
        iss.clear();
        iss.str(*iter++);
        iss >> data.numFrags_;
        lineOK = lineOK && !iss.fail();
    }

    // Read remaining fragments
    data.frags_.clear();
    data.frags_.reserve(data.numFrags_);
    while(iter != E && lineOK)
    {
        int size;
        iss.clear();
        iss.str(*iter++);
        iss >> size;
        data.frags_.push_back(FragData(size));
        lineOK = lineOK && !iss.fail();
    }

    lineOK = data.frags_.size() == (size_t) data.numFrags_;

    if (!lineOK)
    {
        data.reset();
        ostringstream msg;
        msg << "Error: Map file formatted incorrectly: " << fileName_ << ".";
        throw(Exception(msg.str()));
        return;
    }

    #if MAPREADER_DEBUG > 0
    std::cout << "Read map: " << data << std::endl;
    #endif
}
