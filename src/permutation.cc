#include <algorithm>
#include <cassert>


#include "permutation.h"
#include "match.h"
#include "dp.h"
#include "mapTypes.h"

using namespace std;

// Declarations
ostream& operator<<(ostream& os, const FloatVec& fv);
istream& operator>>(istream& is, FloatVec& fv);

// Compute the fraction of scores that are as good or better
// than this score for this number of fragments
float PermutationDB::getPval(unsigned int numFrags, float score) const {
    
    UInt2FloatVec::const_iterator p = scoreMap_.find(numFrags);

    if (p == scoreMap_.end()) {
        assert(false);
        return 1.0;
    }

    const FloatVec& fv = p->second;
    assert(fv.size() == numTrials_);
    const FloatVec::const_iterator E = fv.end();
    const FloatVec::const_iterator pos = lower_bound(fv.begin(), E, score);
    int numBetter = E - pos;
    return ((float) numBetter)/numTrials_;
}


// Run Permutation tests for random restriction patters up to
// maxN fragments in length
void PermutationDB::runTests(unsigned int maxN, const AlignmentParams& alignParams, const vector<OpticalMapData *>& opMaps)
{
    scoreMap_.clear();
    // Set default entries in scoreMap_;
    FloatVec infVec = FloatVec(numTrials_, -Constants::INF);
    for(unsigned int i = 1; i <= maxN; i++)
        scoreMap_.insert(UInt2FloatVec::value_type(i, infVec));


    // For each fragment number, perform a permutation test
    UInt2FloatVec::iterator i = scoreMap_.begin();
    UInt2FloatVec::iterator e = scoreMap_.end();
    size_t numRounds = scoreMap_.size();
    size_t numMaps = opMaps.size();
    uint reportPeriod = ((uint) 0.01*numRounds) + 1;
    uint numProcessed = 0;
    for(; i != e; i++)
    {
        if (numProcessed%reportPeriod == 0)
            cout << "Processed permutation test " << numProcessed << " of " << numRounds << endl;

        // Generate random contigs
        uint numFrags = i->first;
        FloatVec& scores = i->second;
        vector< vector<FragData *> > permutedContigs(numTrials_);
        for(uint j = 0; j < numTrials_; j++)
            permutedContigs[j] = generateRandomContig(numFrags);

        #pragma omp parallel for schedule(dynamic) private(j)
        for (uint j=0; j < numTrials_; j++)
        {
            // Match permuted contig to all optical maps
            vector<FragData *>& randomFrags = permutedContigs[j];
            double maxScore = -Constants::INF;
            for (uint k=0; k<numMaps; k++)
            {
                const OpticalMapData * pOpMap = opMaps[k];
                double scoreForward =  matchPermutationTest(randomFrags, pOpMap, alignParams);
                reverse(randomFrags.begin(), randomFrags.end());
                double scoreReverse =  matchPermutationTest(randomFrags, pOpMap, alignParams);
                if (scoreForward > maxScore) maxScore = scoreForward;
                if (scoreReverse > maxScore) maxScore = scoreReverse;
            }
            scores[j] = maxScore;
        }
        sort(scores.begin(), scores.end()); // sort the scores
        numProcessed++;
    }
}

// Run Permutation tests for random restriction patters up to
// maxN fragments in length
void PermutationDB::runTests(const set<unsigned int>& numFrags, const AlignmentParams& alignParams, const vector<OpticalMapData *>& opMaps)
{
    scoreMap_.clear();
    // Set default entries in scoreMap_;
    FloatVec infVec = FloatVec(numTrials_, -Constants::INF);
    {
    set<unsigned int>::const_iterator iter = numFrags.begin();
    set<unsigned int>::const_iterator E = numFrags.end();
    for(; iter != E; iter++)
        scoreMap_.insert(UInt2FloatVec::value_type(*iter, infVec));
    }

    // For each fragment number, perform a permutation test
    UInt2FloatVec::iterator i = scoreMap_.begin();
    UInt2FloatVec::iterator e = scoreMap_.end();
    size_t numRounds = scoreMap_.size();
    size_t numMaps = opMaps.size();
    uint reportPeriod = ((uint) 0.01*numRounds) + 1;
    uint numProcessed = 0;
    for(; i != e; i++)
    {
        if (numProcessed%reportPeriod == 0)
            cout << "Processed permutation test " << numProcessed << " of " << numRounds << endl;

        // Generate random contigs
        uint numFrags = i->first;
        FloatVec& scores = i->second;
        vector< vector<FragData *> > permutedContigs(numTrials_);
        for(uint j = 0; j < numTrials_; j++)
            permutedContigs[j] = generateRandomContig(numFrags);

        #pragma omp parallel for schedule(dynamic) private(j)
        for (uint j=0; j < numTrials_; j++)
        {
            // Match permuted contig to all optical maps
            vector<FragData *>& randomFrags = permutedContigs[j];
            double maxScore = -Constants::INF;
            for (uint k=0; k<numMaps; k++)
            {
                const OpticalMapData * pOpMap = opMaps[k];
                double scoreForward =  matchPermutationTest(randomFrags, pOpMap, alignParams);
                reverse(randomFrags.begin(), randomFrags.end());
                double scoreReverse =  matchPermutationTest(randomFrags, pOpMap, alignParams);
                if (scoreForward > maxScore) maxScore = scoreForward;
                if (scoreReverse > maxScore) maxScore = scoreReverse;
            }
            scores[j] = maxScore;
        }
        sort(scores.begin(), scores.end()); // sort the scores
        numProcessed++;
    }
}
// Populate vector of contig fragments. These are sampled
// to form random contig frags.
void PermutationDB::setContigFrags(const vector<ContigMapData *>& contigVec)
{
    contigFrags_.clear();

    // Collect all contig frags
    vector<ContigMapData *>::const_iterator i = contigVec.begin();
    const vector<ContigMapData *>::const_iterator e = contigVec.end();
    vector<FragData *>::const_iterator fragIter, fragEnd;
    for(; i != e; i++)
    {
       const ContigMapData * pMap = *i;
       vector<FragData *>::const_iterator fragEnd = pMap->frags_.end();
       vector<FragData *>::const_iterator fragStart = pMap->frags_.begin();
       contigFrags_.insert(contigFrags_.end(), fragStart, fragEnd);
    }

    N_ = contigFrags_.size();
}

// Generate a random contig by sampling with replacement from the distribution of 
// contig fragment sizes
vector<FragData *> PermutationDB::generateRandomContig(unsigned int numFrags) const
{
    vector<FragData *> frags(numFrags);
    for(unsigned int k=0; k < numFrags; k++)
        frags[k] =  contigFrags_[rand() % N_]; // random frag
    return frags;
}

void PermutationDB::write(ostream& os) const {

    UInt2FloatVec::const_iterator i = scoreMap_.begin();
    const UInt2FloatVec::const_iterator E = scoreMap_.end();
    for(; i != E; i++) {
       unsigned int numFrags = i->first;
       const FloatVec& fv = i->second;
       os << numFrags << " " << fv << "\n";
    }
}

void PermutationDB::read(istream& is) {
    reset();
    string line;
    uint fragCount;
    FloatVec scores;
    while(getline(is, line))
    {
        istringstream iss(line);
        iss >> fragCount >> scores;
        assert(!is.fail());
        assert(scoreMap_.count(fragCount)==0);
        scoreMap_.insert(UInt2FloatVec::value_type(fragCount, scores));
    }
}

ostream& operator<<(ostream& os, const FloatVec& fv) {
    FloatVec::const_iterator i = fv.begin();
    const FloatVec::const_iterator E = fv.end();
    const FloatVec::const_iterator last = E-1;
    for(; i != E; i++) {
        os << *i;
        if (i != last) os << " ";
    }
    return os;
}

istream& operator>>(istream& is, FloatVec& fv)
{
    fv.clear();
    float num;
    while (is >> num) {
        fv.push_back(num);
    }
    return is;
}
