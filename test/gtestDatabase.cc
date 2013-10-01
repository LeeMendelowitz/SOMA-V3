// Test the FragDatabase
#include "gtest/gtest.h"

#include <iostream>

#include "ContigMapData.h"
#include "OpticalMapData.h"
#include "FragDatabase.h"

// To use a test fixture, derive a class from testing::Test.
class DatabaseTest : public testing::Test {
 protected:  // You should make the members protected s.t. they can be
             // accessed from sub-classes.

  // virtual void SetUp() will be called before each test is run.  You
  // should define it if you need to initialize the varaibles.
  // Otherwise, this can be skipped.
  virtual void SetUp()
  {
    readOpticalMaps("map100.opt", opticalMaps_);
    readContigMaps("contigs.silico", contigMaps_);
    std::cout << "Read " << contigMaps_.size() << " contig maps.\n";
    std::cout << "Read " << opticalMaps_.size() << " optical maps.\n";
    fragDB_.addMap(opticalMaps_.front());
    fragDB_.sortFrags();
  }

  // virtual void TearDown() will be called after each test is run.
  // You should define it if there is cleanup work to do.  Otherwise,
  // you don't have to provide it.
  //
  virtual void TearDown() {
    for(size_t i = 0; i < opticalMaps_.size(); i++)
        delete opticalMaps_[i];
    opticalMaps_.clear();
    for(size_t i = 0; i < contigMaps_.size(); i++)
        delete contigMaps_[i];
    contigMaps_.clear();
  }

 
  void testFragConstructor()
  { }

  // Declares the variables your tests want to use.
  vector<ContigMapData*> contigMaps_;
  vector<OpticalMapData*> opticalMaps_;
  FragDatabase fragDB_;
};


bool hitWithinBounds(const IntPairVec& bounds, FragPtr * hit)
{
    FragDataVec::const_iterator d = hit->pFrag_;
    size_t i = 0;
    bool isMatch = true;
    for(; i < bounds.size(); i++, d++)
    {
        int lb = bounds[i].first;
        int ub = bounds[i].second;
        isMatch = isMatch && (d->size_ <= ub) && (d->size_ >= lb);
        if (!isMatch) break;
    }
    return isMatch;
}

bool indexIsCorrect(FragPtrVec::const_iterator B, FragPtrVec::const_iterator E, int ind, int q)
{

    bool isOK = true;
    FragPtrVec::const_iterator iter = B + ind;

    // Check that this is the first element greater than q
    if (iter < E)
    {
        FragPtr * p = *iter;
        isOK = isOK &&  (p->pFrag_->size_ > q);
    }

    if (iter - 1 >= B)
    {
        FragPtr * p = *(iter-1);
        isOK = isOK && (p->pFrag_->size_ < q);
    }
    return isOK;
}


// When you have a test fixture, you define a test using TEST_F
// instead of TEST.

// Tests the default c'tor.
TEST_F(DatabaseTest, ConstructerSortTest) {
    FragPtrVec::const_iterator i = fragDB_.fragsB();
    int lastSize = (*i)->pFrag_->size_;
    int lastRank = (*i)->rank_;
    ASSERT_EQ(0, lastRank);
    i++;
    for(; i != fragDB_.fragsE(); i++)
    {
       FragPtr * pFrag = *i;
       ASSERT_GE(pFrag->pFrag_->size_, lastSize);
       ASSERT_EQ(lastRank + 1, pFrag->rank_);
       lastSize = pFrag->pFrag_->size_;
       lastRank++;
    }
}

TEST_F(DatabaseTest, FragPtrPrevNextTest) {
    FragPtrVec::const_iterator i = fragDB_.fragsB();
    for(; i != fragDB_.fragsE(); i++)
    {
       FragPtr * pFrag = *i;

       // Test that the next pointer is correct
       FragDataVec::const_iterator next = pFrag->pFrag_ + 1;
       FragDataVec::const_iterator prev = pFrag->pFrag_ - 1;
       if (next < pFrag->map_->getFragsE())
       {
           ASSERT_EQ(next, pFrag->pNext_->pFrag_);
       }
       else
       {
           ASSERT_EQ(NULL, pFrag->pNext_);
       }

       if (prev >= pFrag->map_->getFragsB())
       {
            ASSERT_EQ(prev, pFrag->pPrev_->pFrag_);
       }
       else
       {
           ASSERT_EQ(NULL, pFrag->pPrev_);
       }
    }
}

TEST_F(DatabaseTest, LowerBoundTest)
{
    std::vector<int> ints;
    ints.push_back(2000);
    ints.push_back(5000);
    ints.push_back(10000);
    ints.push_back(30000);
    ints.push_back(50000);
    ints.push_back(10000000);

    for (size_t i = 0; i < ints.size(); i++)
    {
        int q = ints[i];
        int ind = fragDB_.lowerBound(q);
        bool isOK = indexIsCorrect(fragDB_.fragsB(), fragDB_.fragsE(), ind, q);
        ASSERT_TRUE(isOK);
    }

}

TEST_F(DatabaseTest, QueryTest)
{
    // Make a query
    IntPairVec bounds;
    // 40041   33760   49438
    bounds.push_back(IntPair(30000,50000));
    bounds.push_back(IntPair(30000,50000));
    bounds.push_back(IntPair(30000,50000));

    vector<FragPtr*> hits;
    fragDB_.getFragPtrHits(bounds, hits);

    cout << "Found " << hits.size() << " hits.\n";
    for(size_t i = 0; i < hits.size(); i++)
    {
        cout << *(hits[i]) << "\n";
    }

    for(size_t i=0; i < hits.size(); i++)
    {
        ASSERT_TRUE(hitWithinBounds(bounds, hits[i]));
    }
}
