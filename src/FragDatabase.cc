#include <algorithm>
using std::sort;

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
