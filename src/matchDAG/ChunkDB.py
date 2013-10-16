"""
Code to find alignments of restriction maps.
"""
import numpy as np
import sys

from DictWrap import DictWrap
import SOMAMapUtils
from Chunk import chunksFromMap

def debugMsg(msg):
    sys.stderr.write(msg)
    sys.stderr.flush()


class ChunkDB(object):
    """
    Class to efficiently query fragment matches.
    """
    def __init__(self, chunks=[]):
        self.chunks = []
        if chunks:
            self.addChunks(chunks)
            self.sortChunks()

    def addChunks(self, chunks):
        chunks = [c for c in chunks]
        self.chunks.extend(chunks)

    def sortChunks(self):
        debugMsg('ChunkDB: sorting chunks...')
        self.chunks = np.sort(self.chunks)
        debugMsg('done!\n')

    def getChunks(self, lowerBound, upperBound):
        lbi = np.searchsorted(self.chunks, lowerBound, side='left')
        ubi = np.searchsorted(self.chunks, upperBound, side='right')
        return self.chunks[lbi:ubi]

 
def processMapFile(mapFile, maxMisses = 5):
    mapd = DictWrap(SOMAMapUtils.readMaps(mapFile))
    allChunks = [c for m in mapd.itervalues() for c in chunksFromMap(m, maxMisses)]
    chunkDB = ChunkDB(chunks=allChunks)
    ret = { 'mapd' : mapd,
            'chunkDB' : chunkDB}
    return DictWrap(ret)
