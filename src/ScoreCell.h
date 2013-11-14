#ifndef SCORECELL_H
#define SCORECELL_H

// ScoreCell is a cell in the dynamic programming table and
// identifies its successors/predecessors.
// It acts as a node in the DAG underlying the dynamic programming
// table.

#include<vector>
#include<set>
#include <unordered_map>
#include <ostream>

#include "ScorePathStep.h"

class ScoreCell;
typedef std::vector<ScoreCell *> ScoreCellPVec;
typedef std::set<ScoreCell *> ScoreCellSet;
typedef std::pair<int, int> IntPair;

class ScoreCell
{
    public:
    ScoreCell() : q_(-1), r_(-1), inPlay_(false), backPointer_(nullptr, nullptr) { };
    ScoreCell(int q, int r, bool inPlay = false) :
        q_(q), r_(r), inPlay_(inPlay), backPointer_(nullptr, nullptr) { };
    ScoreCell(IntPair ip, bool inPlay = false) :
        q_(ip.first), r_(ip.second), inPlay_(inPlay), backPointer_(nullptr, nullptr) { };

    IntPair key() const { return IntPair(q_, r_);}
    bool operator==(const ScoreCell& other) const { return (key() == other.key()); } 
    size_t inDegree() const { return backPointers_.size(); };
    size_t outDegree() const { return forwardPointers_.size(); };
    bool removeBackPointer(const ScoreCell* dest);
    bool removeForwardPointer(const ScoreCell* dest);
    void removeBackPointers() { backPointers_.clear(); }
    void removeForwardPointers() { forwardPointers_.clear(); }
    void reset();

    int q_;
    int r_;
    bool inPlay_;
    double score_;

    ScorePathStepVec backPointers_;
    ScorePathStepVec forwardPointers_;
    ScorePathStep backPointer_; // back pointer for DP solution path
};

/*
namespace std {
    template <>
    struct std::hash< ScoreCell > {
        public:
        size_t operator()(const ScoreCell& x) const throw() {
            std::hash<IntPair> hasher;
            return hasher(x.key());
        }
    };
}
*/

inline void ScoreCell::reset()
{
    q_ = -1;
    r_ = -1;
    inPlay_ = false;
    backPointers_.clear();
    forwardPointers_.clear();
    backPointer_ = ScorePathStep(nullptr, nullptr);
}

typedef std::vector<ScoreCell> ScoreCellVec; typedef unordered_map<IntPair, ScoreCell *> ScoreCellMap;

std::ostream& operator<<(std::ostream&, const ScoreCell& cell);

#endif
