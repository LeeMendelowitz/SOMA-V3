import sys
from StringIO import StringIO
import numpy as np
from operator import attrgetter
from collections import defaultdict
from datetime import datetime

from Chunk import chunksFromMap

def DEBUG(msg):
    sys.stderr.write(msg)
    sys.stderr.flush()


#################################################################
class ChunkPairNode(object):
    """
    Node in the match search tree.
    """
    def __init__(self, parent, queryChunk, refChunk, numQueryMisses, numRefMisses):
        self.queryChunk = queryChunk
        self.refChunk = refChunk
        self.parent = parent
        self.numQueryMisses = numQueryMisses
        self.numRefMisses = numRefMisses


# Check if chunks are compatible
def chunksCompatible(chunk1, chunk2, minDist, maxDist, minSites, maxSites):
    """
    Two chunks are compatible if they are from the same map, if the distance
    separating them is within the bounds, and if the number of interior sites is within the bounds.
    """
    if not (minDist <= maxDist):
        raise RuntimeError('minDist must be less than maxDist')

    if not (minSites <= maxSites):
        raise RuntimeError('minSites must be less than maxSites')

    if chunk1.map.mapId != chunk2.map.mapId:
        return False

    if chunk1.bInd < chunk2.bInd:
        leftChunk, rightChunk = chunk1, chunk2
    else:
        leftChunk, rightChunk = chunk2, chunk1
    
    dist = rightChunk.leftPos - leftChunk.rightPos
    numSites = rightChunk.bInd - leftChunk.eInd

    return (minDist <= dist <= maxDist) and (minSites <= numSites <= maxSites)

class GraphBuilderDoubleSided(object):

    @staticmethod
    def markRight(chunkPairs):
        for q,r in chunkPairs:
            r._rmark = True
            if hasattr(r, '_qChunks'):
                r._qChunks.append(q)
            else:
                r._qChunks = [q]

    @staticmethod
    def unmarkRight(chunkPairs):
        refChunks = (r for q,r in chunkPairs)
        for chunk in refChunks:
            if hasattr(chunk, '_rmark'):
                del chunk._rmark
            if hasattr(chunk, '_qChunks'):
                del chunk._qChunks

    @staticmethod
    def unmarkRight2(chunkPairs):
        refChunks = (r for q,r in chunkPairs)
        for chunk in refChunks:
            try:
                del chunk_rmark
            except AttributeError:
                pass
            try:
                del chunk._qChunks
            except AttributeError:
                pass

    def getChunkPairs(self, minDist, maxDist, minSites, maxSites):

        # Mark the right reference chunks
        self.markRight(self.endPairs)

        # For each left reference chunk, search for a marked right reference chunk
        # within the distance bounds.


        # Unmark the right reference chunks
        self.unmarkRight(self.endPairs)

    def __init__(self, queryMap, refChunkDB, tol=0.10, minDelta=1000, maxNodes=100000, maxMissRate = 0.5, maxQueryMisses=3, maxRefMisses=3):

        # Generate all chunks from the queryMap
        qChunks = chunksFromMap(queryMap, maxQueryMisses)
        qsChunks = [c for c in qChunks if c.bInd == 0]
        qeChunks = [c for c in qChunks if c.eInd == queryMap.numFrags]
        DEBUG('Have %i query start chunks.\n'%len(qsChunks))
        DEBUG('Have %i query end chunks.\n'%len(qeChunks))

        compatiblePairs = []

        def getRefMatches(q):
            delta = max(minDelta, tol*q.size)
            lb, ub = q.size - delta, q.size + delta
            refChunks = refChunkDB.getChunks(lb, ub)
            return refChunks

        # Compute starting/ending chunks pairs
        DEBUG('Buildling startPairs/endPairs\n')
        startPairs = [(q, r) for q in qsChunks for r in getRefMatches(q)]
        DEBUG('Made %i startPairs\n'%len(startPairs))
        endPairs = [(q,r) for q in qeChunks for r in getRefMatches(q)]
        DEBUG('Made %i endPairs\n'%len(endPairs))
        return

        self.startPairs = startPairs
        self.endPairs = endPairs

        self.qToRStart = defaultdict(list)
        self.qToREnd = defaultdict(list)

        for q,r in self.endPairs:
            if not hasattr(r, '_qEndList'):
                r._qEndList = [q]
            else:
                r._qEndList.append(q)


        # Bin by reference map
        refMapToStartPairs = defaultdict(list)
        refMapToEndPairs = defaultdict(list)
        for qs,rs in startPairs:
            refMapToStartPairs[rs.map].append((qs,rs))
        for qe,re in endPairs:
            refMapToEndPairs[re.map].append((qe,re))
        self.refMapToStartPairs = refMapToStartPairs
        self.refMapToEndPairs = refMapToEndPairs

        allMaps = list(set(refMapToStartPairs.keys() + refMapToEndPairs.keys()))
        for m in allMaps:
            mapStartPairs = refMapToStartPairs[m]
            mapEndPairs = refMapToEndPairs[m]
            DEBUG("map:%s startPairs:%i endPairs:%i\n"%(m.mapId, len(mapStartPairs), len(mapEndPairs)))
            for startPair in mapStartPairs:
                qs,rs = startPair
                for endPair in mapEndPairs:
                    qe,re = endPair
                    numQuerySites = qe.bInd - qs.eInd
                    queryDist =  qe.leftPos - qs.rightPos
                    delta = max(minDelta, tol*queryDist)
                    minRefDist = queryDist - delta
                    maxRefDist = queryDist + delta
                    deltaSites = maxMissRate*numQuerySites
                    minRefSites = numQuerySites - deltaSites
                    maxRefSites = numQuerySites + deltaSites
                    if chunksCompatible(rs, re, minRefDist, maxRefDist, minRefSites, maxRefSites):
                        compatiblePairs.append((startPair, endPair))

        self.compatiblePairs = compatiblePairs
        DEBUG('Found %i compatible pairs.\n'%len(compatiblePairs))


class GraphBuilder(object):

    def __init__(self, queryChunk, refChunkDB, tol=0.10, minDelta=1000, maxNodes=100000, maxMissRate = 0.5):

        start = datetime.now()

        self.roots = []
        self.graphs = []

        # Determine how many missed sites are permitted interior to global alignment
        numQueryFrags = queryChunk.map.numFrags
        maxMisses = maxMissRate*(numQueryFrags + 1)

        # Find matching reference fragments. Each match gives a root
        # for a search tree
        delta = max(minDelta, tol*queryChunk.size)
        lb, ub = queryChunk.size - delta, queryChunk.size + delta
        refChunks = refChunkDB.getChunks(lb, ub)
        DEBUG('Have %i ref chunks\n'%len(refChunks))
        for refChunk in refChunks:
            numQueryMisses = queryChunk.numFrags - 1
            numRefMisses = refChunk.numFrags - 1
            self.roots.append(ChunkPairNode(None, queryChunk, refChunk, numQueryMisses, numRefMisses))

        # For each match, build a search tree.
        for root in self.roots:
            G = Graph(root, maxMisses=maxMisses)
            self.graphs.append(G)
        self.graphs = sorted(self.graphs, key = attrgetter('numNodes'), reverse=True)
        numPaths = sum(g.numLeaves for g in self.graphs)
        DEBUG('Made %i paths\n'%numPaths)

        end = datetime.now()
        delta = (end-start).total_seconds()
        DEBUG("%.2f seconds elapsed\n'"%delta)

    def dumpDebug(self, handle=None):

        opened = False
        if handle is None:
            opened = True
            handle = open('debug.out', 'w')
        elif isinstance(handle, str):
            opened=True
            handle = open(handle, 'w')

        # Dump debug information for each graph
        numPaths = sum(g.numLeaves for g in self.graphs)
        handle.write('numGraphs:%i\t'%len(self.graphs))
        handle.write('numPaths:%i\t'%numPaths)
        handle.write('\n')
        for g in self.graphs:
            g.dumpDebug(handle)
        if opened:
            handle.close()

class Graph(object):

    def __init__(self, root, tol=0.10, minDelta=1000, maxNodes = 100000, maxMisses=0):
        self.root = root
        self.finished = []
        self.dead = []
        self.nodeCount = 1
        self.leafCount = 1
        self.buildFromRoot(root, tol, minDelta, maxNodes, maxMisses)

    @property
    def numLeaves(self):
        return len(self.finished)

    @property
    def numNodes(self):
        return self.nodeCount

    @property
    def numDead(self):
        return len(self.dead)

    def buildFromRoot(self, root, tol, minDelta, maxNodes, maxMisses):
        leaves = [root]
        finished = []
        dead = []
        nodeCount = 1

        depth = 0
        self.depthToNodeCount = {}
        self.depthToNodeCount[0] = 1
        while (nodeCount < maxNodes) and leaves:
            depth += 1
            newLeaves = []
            for leaf in leaves:
                hasChild = False
                for qc in leaf.queryChunk.successors:
                    numInteriorSitesQuery = qc.numInteriorSites
                    delta = max(minDelta, tol*qc.size)
                    lb, ub = qc.size - delta, qc.size + delta
                    refChunkSuccessors = (rc for rc in root.refChunk.successors if lb <= rc.size <= ub)
                    for rc in refChunkSuccessors:
                        numRefMisses = leaf.numRefMisses + rc.numInteriorSites
                        numQueryMisses = leaf.numQueryMisses + numInteriorSitesQuery
                        if (numRefMisses > maxMisses) or (numQueryMisses > maxMisses):
                            continue
                        child = ChunkPairNode(leaf, qc, rc, numQueryMisses, numRefMisses)
                        hasChild = True
                        newLeaves.append(child)

                # If no children were made, store this as a leaf.
                # The leaf is dead if it is not the end of the query.
                # The leaf is finished if we made it to the end of the query.
                if not hasChild:
                    if leaf.queryChunk.successors:
                        dead.append(leaf)
                    else:
                        finished.append(leaf)
            leaves = newLeaves 
            nodeCount += len(leaves)
            if leaves:
                self.depthToNodeCount[depth] = len(leaves)

        if nodeCount >= maxNodes:
            DEBUG("Reached Node Limit!\n")

        self.finished = finished
        self.dead = dead
        self.nodeCount = nodeCount
        self.leafCount = len(self.finished)

    def _buildFromLeaf(self, leaf):
        # Build the return path from a leaf
        path = [leaf]
        curNode = leaf
        while curNode:
            path.append(curNode)
            curNode = curNode.parent
        path = path[::-1]
        return path

    def iterPaths(self):
        for leaf in self.finished:
            yield self._buildFromLeaf(leaf)

    def dumpDebug(self, handle):
        depths = sorted(self.depthToNodeCount.keys())
        countData = ((d,self.depthToNodeCount[d]) for d in depths)
        handle.write('numNodes:%i\t'%self.nodeCount)
        handle.write('leafCount:%i\t'%self.leafCount)
        handle.write('deadLeafCount:%i\t'%len(self.dead))
        handle.write('nodeCountsByDepth:'+'\t'.join('{0};{1}'.format(d,c) for d,c in countData) + '\n')
        handle.write('\n')

class Score(object):
    def __init__(self, refMissScore, queryMissScore, chi2):
        self.refMissScore = refMissScore
        self.queryMissScore = queryMissScore
        self.chi2 = chi2
        self.total = self.refMissScore + self.queryMissScore + self.chi2

    def __lt__(self, other):
        return self.total < other.total

    def __eq__(self, other):
        return self.total == other.total

    def __str__(self):
        s = StringIO()
        w = s.write
        w('Total: %.2f\n'%self.total)
        w('refMissScore: %.2f\n'%self.refMissScore)
        w('queryMissScore: %.2f\n'%self.queryMissScore)
        w('chi2: %.2f\n'%self.chi2)
        out = s.getvalue()
        s.close()
        return out

    def __repr__(self):
        return str(self)
        

class Scorer(object):
    """
    Class to score an alignment path.

    Sizing error for Valouev et al, "Alignment of Optical Maps", 2006:
        X ~ N(y, sigma^2*y)
        where X is the measured fragment size
        y is the true fragment size
        sigma^2 is a parameter controlling the variance of the measurement.
    """
    # Note: stddev(y) = sigma*sqrt(y)
    # The default values below imply stddev(100kb) = 10 kb, stddv(10kb) = 3.1kb

    # THESE DEFAULT VALUES ARE VERY LOOSE. PROBABLY SHOULD CHANGE THEM.
    default_sigma_kb = 1.0
    default_sigma_kb = 0.5
    default_sigma_bp = np.sqrt(1000.0)*default_sigma_kb

    def __init__(self, refMissPenalty = 5, queryMissPenalty = 3, sigma = default_sigma_bp):
        self.refMissPenalty = refMissPenalty
        self.queryMissPenalty = queryMissPenalty
        self.sigma = sigma


    def score(self, path):
        refChunks = [n.refChunk for n in path]
        queryChunks = [n.queryChunk for n in path]
        numRefMisses = sum(c.refChunk.numInteriorSites for c in path)
        numQueryMisses = sum(c.queryChunk.numInteriorSites for c in path)
        rl = np.array([c.size for c in refChunks])
        ql = np.array([c.size for c in queryChunks])
        deltas = ql - rl
        deltas2 = deltas**2
        var = self.sigma**2*rl
        chi2 = np.sum(deltas2/var)

        refMissScore = -numRefMisses*self.refMissPenalty
        queryMissScore = -numQueryMisses*self.queryMissPenalty
        chi2 = -chi2
        score = Score(refMissScore, queryMissScore, chi2)
        score.rl = rl
        score.ql = ql
        score.deltas = deltas
        score.var = var
        score.std = np.sqrt(var)
        score.chi2_comps = deltas2/var
        score.numRefMisses = numRefMisses
        score.numQueryMisses = numQueryMisses
        return score

########################################################################################



class Node3(object):

    def __init__(self, queryChunk, refChunks, parent = None):
        self.queryChunk = queryChunk
        self.refChunks = refChunks
        self.parent = parent


def leftRightIntersect(leftNode, rightNodes):
    """
    Take the intersection by keeping chunks in the left chunkset with
    a successor in the rightChunkSet.
    Return the new left set.
    """
    rightNodeCoordUnion = set()
    for rn in rightNodes:
        rightNodeCoordUnion = rightNodeCoordUnion.union(rn._refChunkCoordSet)

    leftChunks = [c for c in leftNode._refChunks if (c.map, c.eInd) in rightNodeCoordUnion]
    leftNode._refChunks = leftChunks
    leftNode._refChunkCoordSet = set((c.map, c.bInd) for c in leftChunks)
    return leftChunks

class GraphBuilder3(object):

    def addRefChunks(self):
        """
        Add compatible reference chunks to each query chunk.
        """
        queryChunks = self.qChunks
        qToRefMatches = self.qToRefMatches

        for q in queryChunks:
            refChunks = qToRefMatches[q]
            q._refChunks = refChunks
            #q._refChunkCoordSet = set((c.map, c.bInd) for c in refChunks)
            q._refChunkCoordSet = set()

    def delRefChunks(self):
        """
        Remove compatible reference chunks from each query chunk.
        """
        for q in self.qChunks:
            del q._refChunks

    def __init__(self, queryMap, refChunkDB, tol=0.10, minDelta=1000, maxNodes=100000, maxMissRate = 0.5, maxQueryMisses=3, maxRefMisses=3):

        start = datetime.now()

        # Generate all chunks from the queryMap
        qChunks = chunksFromMap(queryMap, maxQueryMisses)
        DEBUG('Have %i query chunks.\n'%len(qChunks))

        def getRefMatches(q):
            delta = max(minDelta, tol*q.size)
            lb, ub = q.size - delta, q.size + delta
            refChunks = refChunkDB.getChunks(lb, ub)
            return refChunks

        # Index query chunks based on where they start/end
        qStartToChunk = defaultdict(list)
        qEndToChunk = defaultdict(list)
        for q in qChunks:
            qStartToChunk[q.bInd].append(q)
            qEndToChunk[q.eInd].append(q)
        qsChunks = qStartToChunk[0]
        qeChunks = qEndToChunk[queryMap.numFrags]
        self.qChunks = qChunks

        # For each chunk in the query, find all matches in the reference
        DEBUG('Aligning each query chunk...')
        refMatches = [getRefMatches(q) for q in qChunks]
        DEBUG('DONE.\n')

        qToRefMatches = dict((q,rl) for q,rl in zip(qChunks, refMatches))

        # Store in self
        self.qStartToChunk = qStartToChunk
        self.qEndToChunk = qEndToChunk
        self.qsChunks = qsChunks
        self.qeChunks = qeChunks
        self.refMatches = refMatches
        self.qToRefMatches = qToRefMatches

        # Build a set on the alignments
        DEBUG('Adding ref chunks...\n')
        #matchSet = set((q,r) for q,rl in zip(qChunks, refMatches) for r in rl)
        #matchD = dict((q,set(rl)) for q,rl in zip(qChunks, refMatches))
        self.addRefChunks()
        DEBUG('DONE!\n')
        numPlacementsStart = sum(len(q._refChunks) for q in qsChunks)
        DEBUG('Have %i placements for start chunks.\n'%numPlacementsStart)

        # Sort chunks from left to right based on where they end.
        qChunksSorted = sorted(qChunks, key = lambda c: (c.eInd, c.bInd), reverse=True)

        for qChunk in qChunksSorted:
            if not qChunk.successors:
                continue
            qChunk._refChunks = leftRightIntersect(qChunk, qChunk.successors)

        numPlacementsStart = sum(len(q._refChunks) for q in qsChunks)
        DEBUG('After pruning: Have %i placements for start chunks.\n'%numPlacementsStart)

        self.roots = []
        self.graphs = []

        end = datetime.now()
        delta = (end-start).total_seconds()
        DEBUG('%.2f seconds have elapsed\n'%(delta))

        # Build search graph for each 

        self.graphs = sorted(self.graphs, key = attrgetter('numNodes'), reverse=True)
        numPaths = sum(g.numLeaves for g in self.graphs)
        DEBUG('Made %i paths\n'%numPaths)



    def buildGraph(self, queryChunk, refChunkDB, tol=0.10, minDelta=1000, maxNodes=100000, maxMissRate = 0.5):

        start = datetime.now()

        # Determine how many missed sites are permitted interior to global alignment
        numQueryFrags = queryChunk.map.numFrags
        maxMisses = maxMissRate*(numQueryFrags + 1)

        roots = []
        for refChunk in queryChunk._refChunks:
            numQueryMisses = queryChunk.numFrags - 1
            numRefMisses = refChunk.numFrags - 1
            roots.append(ChunkPairNode(None, queryChunk, refChunk, numQueryMisses, numRefMisses))

        self.roots.extend(roots)

        # For each match, build a search tree.
        graphs = []
        for root in roots:
            G = Graph(root, maxMisses=maxMisses)
            graphs.append(G)

        self.graphs.extend(graphs)

        end = datetime.now()
        delta = (end-start).total_seconds()
        DEBUG("%.2f seconds elapsed\n'"%delta)

    def dumpDebug(self, handle=None):

        opened = False
        if handle is None:
            opened = True
            handle = open('debug.out', 'w')
        elif isinstance(handle, str):
            opened=True
            handle = open(handle, 'w')

        # Dump debug information for each graph
        numPaths = sum(g.numLeaves for g in self.graphs)
        handle.write('numGraphs:%i\t'%len(self.graphs))
        handle.write('numPaths:%i\t'%numPaths)
        handle.write('\n')
        for g in self.graphs:
            g.dumpDebug(handle)
        if opened:
            handle.close()

    def getStartPlacements(self):
        return set(c for q in self.qsChunks for c in q._refChunks)

    def buildFromRoot(self, root):
        leafs = [root]
        done = False

        while not done:

            for l in leafs:

                # Get each query successor of the queryChunk
                for qNext in l.queryChunk.successors:

                    rightRefChunkSet = qNext._refPredecessorChunkSet
                    # Take the in




