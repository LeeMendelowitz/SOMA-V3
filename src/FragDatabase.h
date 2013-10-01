/*
The FragDatabase stores reference fragments in a sorted array, which makes it possible
to efficiently find hits for a query.
*/
#ifndef FRAGDATABASE_H
#define FRAGDATABASE_H

#include <vector>

#include "OpticalMapData.h"
#include "seedTypes.h"

typedef std::vector<FragPtr*> FragPtrVec;
typedef std::pair<int, int> IntPair;
typedef std::vector<IntPair> IntPairVec;

class FragDatabase
{

    public:

    FragDatabase() {};
    ~FragDatabase();

    void addMap(const OpticalMapData * pMap);
    void sortFrags();
    int lowerBound(int q);
    bool getFragPtrHits(IntPairVec& query, std::vector<FragPtr*>& hits);

    FragPtrVec::const_iterator fragsB() const { return frags_.begin();}
    FragPtrVec::const_iterator fragsE() const { return frags_.end();}

    private:

    FragPtrVec frags_;
};

#endif
