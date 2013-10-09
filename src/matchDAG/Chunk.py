from IntervalTree import Interval
from collections import defaultdict
from functools import total_ordering
import sys

def makeCmpFuncFromAttrs(attrs):

    def cmpFunc(obj1, obj2):
        for a in attrs:
            v1 = getattr(obj1, a)
            v2 = getattr(obj2, a)
            if v1 < v2:
                return True
            elif v2 < v1:
                return False
        return True
    return cmpFunc

####################################################################################
@total_ordering
class Chunk(object):
    """
    Represents consecutive fragments of a restriction map.
    """
    def __init__(self, map, bInd, eInd, size, isForward=True):
        self.map = map
        self.bInd = bInd # inclusive
        self.eInd = eInd # exclusive
        self.size = size
        self.isForward = isForward # ignored for now
        self.successors = [] # List of Chunks from the same map which immediately follow this one.
        self.predecessors = [] # List of Chunks from the same map which immediately preceed this one.

    @property
    def numFrags(self):
        return self.eInd - self.bInd

    @property
    def numInteriorSites(self):
        return self.numFrags - 1

    def __str__(self):
        fields = [self.map.mapId, self.bInd, self.eInd, self.size]
        fields = [str(f) for f in fields]
        return '(' + ','.join(fields) + ')'

    def __repr__(self):
        return 'Chunk: %s'%(str(self))

    def __lt__(self, other):
        if isinstance(other, (int, float, long)):
            return self.size < other
        elif isinstance(other, Chunk):
            return self.size < other.size
        else:
            raise TypeError('other should be int or Chunk')

    def __eq__(self, other):
        if isinstance(other, (int, float, long)):
            return self.size == other
        elif isinstance(other, Chunk):
            return self.size == other.size
            
    def getSizeOf(self):
        from sys import getsizeof
        attrs = ['map', 'bInd', 'eInd', 'size', 'isForward', 'successors', 'predecessors']
        return sum(getsizeof(getattr(self, a)) for a in attrs)


class ChunkInterval(Interval):
    """
    An interval of sizes which are compatible with the size of the Chunk.
    """
    def __init__(self, chunk, error=0.05, minDelta=1000):   

        self.chunk = chunk

        delta = error*self.chunk.size
        if minDelta:
            delta = max(delta, minDelta)

        # Make interval of width 2*delta
        b = max(0, self.chunk.size - delta)
        e = self.chunk.size + delta
        Interval.__init__(self, b, e)

####################################################################
def _genChunksFromMap(m, maxMisses):
    """
    Generate chunks from a map by considering at most maxMisses missed interior sites.
    """ 
    N = len(m.frags)
    for i in xrange(N):
        s = 0
        for j in xrange(i+1, min(i+maxMisses+2, N)):
            s += m.frags[j-1]
            #s = sum(m.frags[i:j])
            C = Chunk(m, i, j, s)
            yield C

def chunksFromMap(m, maxMisses):
    chunks = [c for c in _genChunksFromMap(m, maxMisses)]
    assignNeighbors(chunks)
    return chunks

##################################################################
def assignNeighbors(chunkList):
    """
    Assign predecessors/successors to Chunks that are in the provided chunkList
    """
    chunkStartD = defaultdict(list)
    chunkEndD = defaultdict(list)
    chunkToStartKey = lambda c: (c.map, c.bInd)
    chunkToEndKey = lambda c: (c.map, c.eInd)

    for c in chunkList:
        chunkStartD[chunkToStartKey(c)].append(c)
        chunkEndD[chunkToEndKey(c)].append(c)
        c.predecessors = []
        c.successors = []

    # In second pass, assign neighbors.
    for c in chunkList:
        chunkEnd = chunkToEndKey(c)
        successors = chunkStartD[chunkEnd]
        c.successors = successors
        for cs in successors:
            cs.predecessors.append(c)

    # Check that everything is sane
    '''
    sys.stderr.write('Checking that everything is sane...')
    sys.stderr.flush()
    for c in chunkList:
       for s in c.successors:
           assert(c in s.predecessors)
       for p in c.predecessors:
           assert(c in p.successors)
    sys.stderr.write('OK!\n')
    sys.stderr.flush()
    '''
