#include <algorithm>
using std::sort;
using std::lower_bound;
using std::upper_bound;
#include <cassert>

#include "FragDatabase.h"

/*
Add fragments from the specified OpticalMapData to the 
database.
*/
void FragDatabase::addMap(const OpticalMapData * pMap)
{

    // Create FragPtr's and put them in the frags_ vector.
    FragPtr * prevFragPtr = 0;
    for (FragDataVec::const_iterator iter = pMap->getFragB();
         iter != pMap->getFragE();
         iter++)
    {

        FragPtr f = FragPtr(pMap, iter);
        frags_.push_back(f);

        FragPtr * curFragPtr = &frags_.back();
        if (prevFragPtr != 0)
        {
            curFragPtr->pPrev_ = prevFragPtr;
            prevFragPtr->pNext_ = curFragPtr;
        }

        prevFragPtr = curFragPtr;
    }
}

FragDatabase::~FragDatabase()
{

    // Delete allocated memory

    // 
    /*
    for (FragPtrVec::iterator iter = frags_.begin();
         iter != frags_.end();
         iter++)
    {

    }
    */
}

///////////////////////////////////////////////////////////////////////////////
/*
    sortFrags

    Sort the fragments, and assign the ranks.

*/
void FragDatabase::sortFrags()
{
    // Sort the fragments
    sort(frags_.begin(), frags_.end());

    // Assign ranks
    int rank = 0;
    for (FragPtrVec::iterator iter = frags_.begin();
         iter != frags_.end();
         iter++)
    {
        iter->rank_ = rank++;
    }
}

inline int FragDatabase::lowerBound(int q)
{
    FragPtrVec::const_iterator iter = lower_bound(frags_.begin(), frags_.end(), q);
    int index = iter - frags_.begin();
    return index;
}

bool FragDatabase::getFragPtrHits(IntPairVec& query, FragPtrVec& hits)
{

    
    for (IntPairVec::const_iterator iter = query.begin();
         iter != query.end();
         iter++)
    {
        int lb = iter->first; // lower bound of query
        int ub = iter->second; // upper bound of query
        assert(lb <= ub);
        int lbi = lowerBound(lb);
        int ubi = lowerBound(ub);
        assert(lbi <= ubi);

        // Use lower index and upper index to get the fragments that match.

        // With each subsequent iteration, move these fragment pointers forward, and filter out those which do not match the next range of indices.

    }


}
