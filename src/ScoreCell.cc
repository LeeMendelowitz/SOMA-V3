#include "ScoreCell.h"

#include <algorithm>



bool ScoreCell::removeBackPointer(const ScoreCell *dest)
{
    /*
    // COMPILE ERROR! Fix if function is to be used
    const ScorePathStepVec::iterator B = backPointers_.begin();
    const ScorePathStepVec::iterator E = backPointers_.end();
    ScorePathStepVec::iterator newEnd = std::remove(B, E, dest);
    bool retVal = (newEnd != E);
    backPointers_.resize(newEnd - B);
    return retVal;
    */
    return false;
}

bool ScoreCell::removeForwardPointer(const ScoreCell *dest)
{
    /*
    // COMPILE ERROR! Fix if function is to be used
    const ScorePathStepVec::iterator B = forwardPointers_.begin();
    const ScorePathStepVec::iterator E = forwardPointers_.end();
    ScorePathStepVec::iterator newEnd = std::remove(B, E, dest);
    bool retVal = (newEnd != E);
    forwardPointers_.resize(newEnd - B);
    return retVal;
    */
    return false;
}

std::ostream& operator<<(std::ostream& os, const ScoreCell& cell)
{
    os << "[(" << cell.q_ << "," << cell.r_ << "), " << cell.score_ << ",";
    const ScoreCell * pTgt = cell.backPointer_.getTarget();
    if (pTgt)
        os << " tgt: (" << pTgt->q_ << "," << pTgt->r_ << ") ]";
    else
        os << " tgt: nullptr ]";
    return os;
}
