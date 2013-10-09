import sys
from StringIO import StringIO
import numpy as np
from operator import attrgetter

from Chunk import chunksFromMap

def DEBUG(msg):
    sys.stderr.write(msg)
    sys.stderr.flush()


###########################################################################################
# OLD CODE

def buildGraph(queryMap, refChunkDB, maxMissesQuery=1, maxDepth=50, maxNodes=1000000):

    # Build chunks on queryMap
    queryChunks = chunksFromMap(queryMap, maxMisses = maxMissesQuery)

    # Get Chunks that start with the beginning of the query.
    startChunks = [c for c in queryChunks if c.bInd == 0]

    startNodes = [Node.makeRoot(c, refChunkDB) for c in startChunks]

    for s in startNodes:
        buildGraphFromRoot(s, maxDepth=maxDepth, maxNodes=maxNodes)

    return startNodes

def bfsSummarize(root):
    print root.summary()
    curNodes = [root]
    depth = 0
    while curNodes:
        print '\n\n\n\n'
        print '*'*50
        print 'Depth: %i Nodes: %i'%(depth, len(curNodes))
        newNodes = []
        for n in curNodes:
            print n.summary()
            newNodes.extend(n.children)
        curNodes = newNodes
        depth += 1
        sys.stdout.flush()
        sys.stdout.write('[ENTER]\n')
        sys.stdout.flush()
        sys.stdin.readline()


def buildGraphFromRoot(root, maxDepth=20, maxNodes=1000000):
    depth = 0    
    nodeCount = 1
    leaves = [root]

    while (depth < maxDepth) and (nodeCount < maxNodes) and leaves:
        DEBUG("*"*50 + '\n')
        DEBUG("Cur depth: %i Num Leaves: %i\n"%(depth, len(leaves)))

        newLeaves = []
        for l in leaves:
            l.makeChildren()
            newLeaves.extend(l.children)
            nodeCount += len(l.children)
            if nodeCount > maxNodes:
                DEBUG("Exceeded max node count! Have %i nodes\n"%nodeCount)
                break
        leaves = newLeaves
        depth += 1

    if depth >= maxDepth:
        DEBUG("Reached max depth %i\n"%depth)


class Node(object):
    """
    Node in the match search tree.
    """
    def __init__(self):
        pass
       
    @staticmethod
    def makeRoot(queryChunk, refChunkDB, tol=0.10, minDelta=1000):
        N = Node()
        N.queryChunk = queryChunk

        # Query the refChunkDB for the compatible reference chunks.
        delta = max(minDelta, tol*queryChunk.size)
        lb, ub = queryChunk.size - delta, queryChunk.size + delta
        N.refChunks = refChunkDB.getChunks(lb, ub)
        N.children = []
        return N

    @staticmethod
    def makeChild(queryChunk, refChunks):
        N = Node()
        N.queryChunk = queryChunk
        N.refChunks = refChunks
        N.children = []
        return N

    def makeChildren(self, tol=0.10, minDelta=1000):
        self.children = []

        """
        DEBUG('Making children for %s\n'%self.queryChunk)
        DEBUG('Query has %i successors\n'%len(self.queryChunk.successors))
        DEBUG('Query has %i reference chunks\n'%len(self.refChunks))
        """

        for querySuccessor in self.queryChunk.successors:
            DEBUG("Working query successor: %s\n"%str(querySuccessor))
            delta = max(minDelta, tol*querySuccessor.size)
            lb = querySuccessor.size - delta
            ub = querySuccessor.size + delta
            compatibleRefs = [refChunkSuccessor for refChunk in self.refChunks for refChunkSuccessor in refChunk.successors
                              if lb <= refChunkSuccessor.size <= ub]
            DEBUG("Found %i compatible refs\n"%len(compatibleRefs))
            if compatibleRefs:
                self.children.append(Node.makeChild(querySuccessor, compatibleRefs))
        DEBUG("Made %i children\n"%len(self.children))
 
    def summary(self):
        s = StringIO()
        w = s.write
        w('Node Summary:\n')
        w('Query Fragment: %s\n'%self.queryChunk)
        w('Num. Ref Chunks: %i\n'%len(self.refChunks))
        w('Num. of Children: %i\n'%len(self.children))
        out = s.getvalue()
        s.close()
        return out

# DONE OLD CODE
##########################################################################

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

class GraphBuilder(object):

    def __init__(self, queryChunk, refChunkDB, tol=0.10, minDelta=1000, maxNodes=100000, maxMissRate = 0.5):

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
