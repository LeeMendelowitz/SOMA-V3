#include <algorithm>
using std::sort;
using std::lower_bound;
using std::upper_bound;
#include <cassert>
#include <iostream>

#include "FragDatabase.h"

void incrementFragPtrs(vector<FragPtr*>& v);
void decrementFragPtrs(vector<FragPtr*>& v);
void filterNullFragPtrs(vector<FragPtr*>& v);
void filterFragPtrsByRank(vector<FragPtr*>&v, int lowerBound, int upperBound);
void filterFragPtrsBySize(vector<FragPtr*>&v, int lowerSize, int upperSize);

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

        FragPtr * f = new FragPtr(pMap, iter);

        if (prevFragPtr != 0)
        {
            f->pPrev_ = prevFragPtr;
            prevFragPtr->pNext_ = f;
        }
        frags_.push_back(f);

        prevFragPtr = f;
    }
}

FragDatabase::~FragDatabase()
{

    // Delete allocated memory
    for (FragPtrVec::iterator iter = frags_.begin();
         iter != frags_.end();
         iter++)
    {
        delete *iter;
    }

    frags_.clear();

}

inline bool FragPtrCmp(const FragPtr * p1, const FragPtr* p2)
{
    return (*p1 < *p2);
}

inline bool FragPtrIntCmp(const FragPtr * p, int q)
{
    return (*p < q);
}

///////////////////////////////////////////////////////////////////////////////
/*
    sortFrags

    Sort the fragments, and assign the ranks.

*/
void FragDatabase::sortFrags()
{
    // Sort the fragments
    sort(frags_.begin(), frags_.end(), FragPtrCmp);

    // Assign ranks
    int rank = 0;
    for (FragPtrVec::iterator iter = frags_.begin();
         iter != frags_.end();
         iter++)
    {
        (*iter)->rank_ = rank++;
    }
}

inline int FragDatabase::lowerBound(int q)
{
    FragPtrVec::const_iterator iter = lower_bound(frags_.begin(), frags_.end(), q, FragPtrIntCmp);
    int index = iter - frags_.begin();
    return index;
}

// Increment all FragPtrs in the vector
void incrementFragPtrs(vector<FragPtr*>& v)
{
    for(vector<FragPtr*>::iterator i = v.begin();
        i != v.end();
        i++)
    {
        FragPtr * old = *i;
        *i = old->pNext_;
    }
    filterNullFragPtrs(v);
}

// Decrement all FragPtrs in the vector
void decrementFragPtrs(vector<FragPtr*>& v)
{
    for(vector<FragPtr*>::iterator i = v.begin();
        i != v.end();
        i++)
    {
        FragPtr * old = *i;
        *i = old->pPrev_;
    }
    filterNullFragPtrs(v);
}

// Remove NULL FragPtrs in the vector
void filterNullFragPtrs(vector<FragPtr*>& v)
{
    const FragPtr* pNull = NULL;
    vector<FragPtr*>::iterator newEnd = remove(v.begin(), v.end(), pNull);
    size_t newSize = newEnd - v.begin();
    v.resize(newSize);
}

// Filter out FragPtr's that do not fall within the rank bounds.
void filterFragPtrsByRank(vector<FragPtr*>&v, int lowerBound, int upperBound)
{
    assert(lowerBound <= upperBound);
    vector<FragPtr*>::iterator iter = v.begin();
    const vector<FragPtr*>::iterator end = v.end();
    vector<FragPtr*>::iterator loc = iter;
    while(iter != end)
    {
        FragPtr* p = *iter;
        if (p->rank_ >= lowerBound && p->rank_ < upperBound)
            *(loc++) = p;
        iter++;
    }

    size_t newSize = loc-v.begin();

    v.resize(newSize);
}

// Filter out FragPtr's that do not fall within the size bounds
void filterFragPtrsBySize(vector<FragPtr*>&v, int lowerSize, int upperSize)
{
    assert(lowerSize <= upperSize);
    vector<FragPtr*>::iterator iter = v.begin();
    const vector<FragPtr*>::iterator end = v.end();
    vector<FragPtr*>::iterator loc = v.begin();
    while(iter != end)
    {
        FragPtr* p = *iter;
        if (p->pFrag_->size_ >= lowerSize && p->pFrag_->size_ < upperSize)
        {
            *loc = p;
            loc++;
        }
        iter++;
    }

    size_t newSize = loc - v.begin();
    v.resize(newSize);
}

bool FragDatabase::getFragPtrHits(IntPairVec& query, vector<FragPtr*>& hits)
{

    if (query.empty())
        return false;

    hits.clear();
    const size_t numFrags = query.size();

    // Populate the hits vector using the first set of bounds.
    IntPairVec::const_iterator qi = query.begin();
    int lb = qi->first; // lower bound of query
    int ub = qi->second; // upper bound of query
    assert(lb <= ub);
    int lbi = lowerBound(lb);
    int ubi = lowerBound(ub);
    assert(lbi <= ubi);
    hits.reserve(ubi - lbi + 1);
    for(int i = lbi; i < ubi; i++)
        hits.push_back(frags_[i]);

    qi++;
    for (; qi != query.end(); qi++)
    {

        // Increment all hits forward
        incrementFragPtrs(hits);
        int lb = qi->first; // lower bound of query
        int ub = qi->second; // upper bound of query
        assert(lb <= ub);

        // Filter out FragPtrs based on rank.
        filterFragPtrsBySize(hits, lb, ub);

        if (hits.empty())
            break;
    }

    // Decrement the FragPtr's by number of fragments in the query.
    const size_t sizeBefore = hits.size();
    for (size_t i = 0; i < numFrags-1; i++)
        decrementFragPtrs(hits);
    const size_t sizeAfter = hits.size();
    assert(sizeBefore == sizeAfter);

    return !hits.empty();
}
