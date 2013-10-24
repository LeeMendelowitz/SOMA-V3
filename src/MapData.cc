#include "MapData.h"
#include "MapChunk.h"

void MapData::indexChunks()
{
    startToChunks_.resize(numFrags());
    endToChunks_.resize(numFrags()+1);

    for (auto c: chunks_)
    {
        startToChunks_[c->getStartIndex()].push_back(c);
        endToChunks_[c->getEndIndex()].push_back(c);
    }
}

void MapData::freeChunks()
{
    for(auto c: chunks_) delete c;
}
