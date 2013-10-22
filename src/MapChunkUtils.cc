#include <vector>
#include <iostream>

#include "MapChunkUtils.h"
#include "MapChunk.h"
#include "MapData.h"
#include "mapTypes.h"

using namespace std;

void setMapChunks(MapData * pMap, size_t maxInteriorMisses)
{

    MapChunkVec chunks;
    const size_t numFrags = pMap->numFrags();
    chunks.reserve(numFrags*(maxInteriorMisses+1));
    MapChunkIndex startToChunks(numFrags);
    MapChunkIndex endToChunks(numFrags + 1);

   
    const FragDataVec::const_iterator fragE = pMap->getFragsE();
    for (FragDataVec::const_iterator iter = pMap->getFragsB();
         iter != fragE;
         iter++)
    {
        for (size_t i = 0; i <= maxInteriorMisses; i++)
        {
            if (iter + i + 1 > fragE) break;
            MapChunk * c = new MapChunk(pMap, iter, iter + i + 1);
            chunks.push_back(c);
            startToChunks[c->getStartIndex()].push_back(c);
            endToChunks[c->getEndIndex()].push_back(c);
        }
    }

    // Set the neighbors of the chunks;
    for(MapChunkVec::const_iterator citer = chunks.begin();
        citer != chunks.end();
        citer++)
    {
        MapChunk * pChunk = *citer;
        size_t offset = pChunk->getEndIndex();
        if (offset < numFrags)
        {
            const MapChunkVec& successors = startToChunks[offset];
            pChunk->next_.reserve(successors.size());
            for(MapChunkVec::const_iterator siter = successors.begin();
                siter != successors.end();
                siter++)
            {
                MapChunk * nextChunk = *siter;
                pChunk->next_.push_back(nextChunk);
                nextChunk->prev_.push_back(pChunk);
            }
        }
    }
    pMap->setChunks(chunks, startToChunks, endToChunks);
}
