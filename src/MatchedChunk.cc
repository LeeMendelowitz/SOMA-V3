#include "MatchedChunk.h"
#include "MapData.h"

// Constructor
MatchedChunk::MatchedChunk(int opStart, int opEnd, int opStartBp, int opEndBp, const MapData* pOpticalMap,
                           int cStart, int cEnd, int cStartBp, int cEndBp, const MapData* pContigMap,
                           bool isBoundaryChunk) :
                           isBoundaryChunk_(isBoundaryChunk)
{
    // Set optical coordinates
    opMapCoord_ = MapCoord(pOpticalMap, opStart, opEnd);
    opSeqCoord_ = SeqCoord(pOpticalMap, opStartBp, opEndBp);

    // Set contig coordinates
    cMapCoord_ = MapCoord(pContigMap, cStart, cEnd);
    cSeqCoord_ = SeqCoord(pContigMap, cStartBp, cEndBp);
}

//
int MatchedChunk::getNumOpticalMisses() const
{
    return opMapCoord_.getNumFrags() - 1;
}

//
int MatchedChunk::getNumContigMisses() const
{
    return cMapCoord_.getNumFrags() - 1;
}

//
int MatchedChunk::getOpticalStartBp() const
{
    return opSeqCoord_.getStartBp();
}

// 
int MatchedChunk::getOpticalEndBp() const
{
    return opSeqCoord_.getEndBp();
}

//
int MatchedChunk::getOpticalMatchLengthBp() const
{
    return opSeqCoord_.getLength();
}
//
int MatchedChunk::getContigStartBp() const
{
    return cSeqCoord_.getStartBp();
}

// 
int MatchedChunk::getContigEndBp() const
{
    return cSeqCoord_.getEndBp();
}

//
int MatchedChunk::getContigMatchLengthBp() const
{
    return cSeqCoord_.getLength();
}

//
int MatchedChunk::getOpticalStartIndex() const
{
    return opMapCoord_.getStart();
}

//
int MatchedChunk::getOpticalEndIndex() const
{

    return opMapCoord_.getEnd();
}

//
int MatchedChunk::getNumOpticalFrags() const
{
    return opMapCoord_.getNumFrags();
}

//
int MatchedChunk::getContigStartIndex() const
{
    return cMapCoord_.getStart();
}

//
int MatchedChunk::getContigEndIndex() const
{

    return cMapCoord_.getEnd();
}

//
int MatchedChunk::getNumContigFrags() const
{
    return cMapCoord_.getNumFrags();
}

//
// Check if this match includes the first or last fragment in a contig
bool MatchedChunk::isBoundaryChunk() const
{
    return isBoundaryChunk_;

//    const MapData * pContigMap = cMapCoord_.getMap();
//    const vector<FragData>& cFrags = pContigMap->getFrags();
//    return (cMapCoord_.getStart() == 0) || ((size_t) cMapCoord_.getEnd() == cFrags.size());
}

//
// Check if this match is an alignment of a contig fragment
// to a gap. 
bool MatchedChunk::isContigGap() const
{
    return (opMapCoord_.getStart() == opMapCoord_.getEnd());
}

//
FragConstIter MatchedChunk::getContigFragB() const
{
    return cMapCoord_.getB();
}

//
FragConstIter MatchedChunk::getContigFragE() const
{
    return cMapCoord_.getE();
}

//
FragConstIter MatchedChunk::getOpticalFragB() const
{
    return opMapCoord_.getB();
}

//
FragConstIter MatchedChunk::getOpticalFragE() const
{
    return opMapCoord_.getE();
}
