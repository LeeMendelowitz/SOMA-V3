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
