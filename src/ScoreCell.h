#ifndef SCORECELL_H
#define SCORECELL_H
// ScoreCell is a cell in the dynamic programming table and
// identifies its successors/predecessors.

#include<vector>

class ScoreCell;
typedef std::vector<ScoreCell *> ScoreCellPVec;

class ScoreCell
{
    public:
    ScoreCell() : q_(-1), r_(-1), inPlay_(false) { };
    ScoreCell(int q, int r, bool inPlay = false) :
        q_(q), r_(r), inPlay_(inPlay) { };
    ScoreCell(IntPair ip, bool inPlay = false) :
        q_(ip.first), r_(ip.second), inPlay_(inPlay) { };

    IntPair key() const { return IntPair(q_, r_);}
    bool operator==(const ScoreCell& other) const { return (key() == other.key()); } 
    size_t inDegree() const { return backPointers_.size(); };
    size_t outDegree() const { return forwardPointers_.size(); };
    void reset();

    int q_;
    int r_;
    bool inPlay_;

    ScorePathStepVec backPointers_;
    ScorePathStepVec forwardPointers_;
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

typedef std::vector<ScoreCell> ScoreCellVec;
typedef unordered_map<IntPair, ScoreCell *> ScoreCellMap;

#endif
