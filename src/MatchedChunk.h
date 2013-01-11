// Represents an aligned fragment block ("chunk")

#ifndef MATCHEDCHUNK_H
#define MATCHEDCHUNK_H

#include "MapCoord.h"
#include "SeqCoord.h"
#include "Scorer.h"

#include <vector>

class MapData;

typedef std::vector<FragData>::const_iterator FragConstIter;

class MatchedChunk
{
    public:

    // Constructor
    MatchedChunk(int opStart, int opEnd, int opStartBp, int opEndBp, const MapData* pOpticalMap,
                int cStart, int cEnd, int cStartBp, int cEndBp, const MapData* pContigMap);

    void setScore(const Score& score) {score_ = score;}

    Score getScore() const { return score_; }
    int getNumOpticalMisses() const;
    int getNumContigMisses() const;
    int getOpticalStartBp() const;
    int getOpticalEndBp() const;
    int getOpticalMatchLengthBp() const;
    int getContigStartBp() const;
    int getContigEndBp() const;
    int getContigMatchLengthBp() const;
    
    int getOpticalStartIndex() const;
    int getOpticalEndIndex() const;
    int getNumOpticalFrags() const;

    int getContigStartIndex() const;
    int getContigEndIndex() const;
    int getNumContigFrags() const;

    bool isBoundaryChunk() const;
    bool isContigGap() const;

    // Iterators to the underlying 
    FragConstIter getContigFragB() const;
    FragConstIter getContigFragE() const;
    FragConstIter getOpticalFragB() const;
    FragConstIter getOpticalFragE() const;

    private:
    MapCoord opMapCoord_; // optical match coordinates
    MapCoord cMapCoord_; // contig match coordinates
    SeqCoord opSeqCoord_; // optical bp coordinates
    SeqCoord cSeqCoord_; // contig bp coordinates
    Score score_; // alignment score
};

#endif
