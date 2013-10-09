from SOMAMapUtils import readMaps
import Chunk
from Match import ChunkDB

refMapFile = 'chr1.map'
queryMapFile = 'MM52_Human_NtBspQI.head.opt'

refMap = readMaps(refMapFile).values()[0]
queryMap = readMaps(queryMapFile).values()[0]

maxRefMisses = 3
maxQueryMisses = 3

refChunks = Chunk.chunksFromMap(refMap, maxRefMisses)
queryChunks = Chunk.chunksFromMap(queryMap, maxQueryMisses)
startChunks = [c for c in queryChunks if c.bInd == 0]
chunkDB =  ChunkDB(refChunks)
queryChunk = startChunks[0]

