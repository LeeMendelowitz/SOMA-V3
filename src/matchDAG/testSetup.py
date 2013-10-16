from SOMAMapUtils import readMaps
import Chunk
from ChunkDB import ChunkDB

refMapFile = 'chr1.map'
refMapFile = 'human.NtBspQI.map'
queryMapFile = 'MM52_Human_NtBspQI.head.opt'

mapD = readMaps(refMapFile)
refMaps = mapD.values()
chr1Map = mapD['chr1.NtBspQI']
queryMap = readMaps(queryMapFile).values()[0]

maxRefMisses = 3
maxQueryMisses = 3
mapChunkList = []
allRefChunks = []
for refMap in refMaps:
    refChunks = Chunk.chunksFromMap(refMap, maxRefMisses)
    mapChunkList.append((refMap, refChunks))
    allRefChunks.extend(refChunks)
mapToChunks = dict((m.mapId, chunks) for m,chunks in mapChunkList)

chr1Chunks = mapToChunks['chr1.NtBspQI']
chr1ChunkDB = ChunkDB(chr1Chunks)

queryChunks = Chunk.chunksFromMap(queryMap, maxQueryMisses)
startChunks = [c for c in queryChunks if c.bInd == 0]
chunkDB =  ChunkDB(allRefChunks)
queryChunk = startChunks[0]

import numpy as np
tol = 0.10
minDelta = 1000
chunkSizes = np.array([q.size for q in queryChunks])
deltas = tol*chunkSizes
deltas = np.array([max(d, minDelta) for d in deltas])
minVals = chunkSizes - deltas
maxVals = chunkSizes + deltas

print 'Querying all chunks'
queryRes = [chunkDB.getChunks(minVal, maxVal) for minVal, maxVal in zip(minVals, maxVals)]
