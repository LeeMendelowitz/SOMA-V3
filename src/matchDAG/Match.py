"""
Perform a global alignment of a query from left to right
using the ChunkDB.

Do not allow any missed sites in the query (i.e. each query site
must be aligned to a reference site.)

This file has some methods for investigating the complexity of the alignment problem.
"""
from DictWrap import DictWrap
import numpy as np
from Clock import Clock
from collections import defaultdict
import sys
from Chunk import chunksFromMap

def DEBUG(msg):
    sys.stderr.write(msg + '\n')
    sys.stderr.flush()



def match(queryMap, chunkDB, minDelta = 1000, tol = 0.1):
    """
    Match the query map to the reference by using the chunkDB
    to align the first fragment in the map, and then iteratively refine.

    Compute the reference chunks where possible alignments end, but do not
    enumerate all possible alignments.

    Ignore the boundary fragments.
    """

    # construct lower bounds and upper bounds
    frags = np.array(queryMap.frags[1:-1])
    numFrags = frags.shape[0]

    # Compute the delta for each fragment
    deltas = tol * frags
    minDeltas = minDelta * np.ones(numFrags)
    deltas = np.max([deltas, minDeltas], axis=0)

    lbs = frags - deltas
    ubs = frags + deltas

    # Query the first fragment
    counts = np.zeros(numFrags)
    expansionCounts = np.zeros(numFrags)
    curRefChunks = chunkDB.getChunks(lbs[0], ubs[0])
    counts[0] = len(curRefChunks)

    for i, (lb,ub) in enumerate(zip(lbs[1:], ubs[1:])):
        #assert(lb <= ub)
        #print 'i=%i'%i
        curRefChunks = set(s for crc in curRefChunks for s in crc.successors)
        expansionCounts[i+1] = len(curRefChunks)
        #print 'After expanding: %i'%len(curRefChunks)
        curRefChunks = [c for c in curRefChunks if lb <= c.size <= ub]
        #print 'After filtering: %i'%len(curRefChunks)
        counts[i+1] = len(curRefChunks)

        if not curRefChunks:
            break

    ret = { 'counts' : counts, 
            'expansionCounts' : expansionCounts,
            'chunks' : curRefChunks }

    return DictWrap(ret)

def dynamicProgrammingTest(queryMap, chunkDB, minDelta=1000, tol=0.1, maxInteriorSitesQuery=3):
    """
    Compute the number of dynamic programming cells that are populated for a given query Map
    """

    # construct lower bounds and upper bounds
    queryChunks = chunksFromMap(queryMap, maxInteriorSitesQuery)

    # Ignore boundary chunks
    numFrags = len(queryMap.frags)
    queryChunks = [c for c in queryChunks if (c.bInd != 0) and (c.eInd != numFrags)]
    chunkSizes = np.array([c.size for c in queryChunks]) 
    numChunks = len(queryChunks)

    # Compute the upper/lower bounds for each fragment
    deltas = tol * chunkSizes
    minDeltas = minDelta * np.ones(numChunks)
    deltas = np.max([deltas, minDeltas], axis=0)
    lbs = chunkSizes - deltas
    ubs = chunkSizes + deltas

    chunkHits = [chunkDB.getChunks(lb, ub) for lb,ub in zip(lbs, ubs)]

    # The coordinates (fragIndex, chromosome name, chromosome chunk end index) give the coordinates
    # of a single cell in the dynamic programming table.
    # For each chromosome, compute the total number of cells, and number of cells that have a match.
    coordGen = ((c.bInd, c.eInd, rc.map, rc.bInd, rc.eInd) for c, refHitList in zip(queryChunks,chunkHits) for rc in refHitList)

    rowToHits = lambda : [set() for i in range(numFrags+1)]
    refMapToHits = defaultdict(rowToHits)
    refMapToInboundHit = defaultdict(rowToHits) # Backpointer points to cell
    refMapToOutboundHit = defaultdict(rowToHits) # Backpointer points from cell
    for qb, qe, refMap, rb, re in coordGen:
        # just add the ending coordinate of chunk
        refMapToHits[refMap][qe].add(re)
        refMapToHits[refMap][qb].add(qe)
        refMapToInboundHit[refMap][qb].add(rb)
        refMapToOutboundHit[refMap][qe].add(re)
        #refMapToHits[refMap].add((qe, re))
        #refMapToHits[refMap].add((qb,qe))

    mapToData = {}
    refMaps = refMapToHits.keys()
    for refMap in refMaps:
        rowToHits = refMapToHits[refMap]
        rowToInboundHit = refMapToInboundHit[refMap]
        rowToOutboundHit = refMapToOutboundHit[refMap]
        numMapFrags = refMap.numFrags
        numQueryFrags = numFrags

        # Count interior cells (i.e. non boundary cases)
        # That have both an inbound and outbound hit.
        maxMisses = maxInteriorSitesQuery
        interiorRows = range(maxMisses+1 , numQueryFrags-maxMisses)
        numInteriorRows = len(interiorRows)
        numBothInboundOutbound = 0
        for irow in interiorRows:
            numBothInboundOutbound += len(rowToInboundHit[irow].intersection(rowToOutboundHit[irow]))


        numCells = numMapFrags*numQueryFrags
        numInteriorCells = numMapFrags*numInteriorRows
        numHits = sum(len(rowHits) for rowHits in rowToHits)
        numInteriorHits = sum(len(rowToHits[r]) for r in interiorRows)
        fracHits = float(numHits)/float(numCells)
        d = {'numMapFrags' : numMapFrags,
             'numQueryFrags' : numQueryFrags,
             'numQueryChunks' : len(queryChunks),
             'numCells' : numCells,
             'numInteriorCells' : numInteriorCells,
             'numInteriorHits' : numInteriorHits,
             'numHits' : numHits,
             'fracHits' : fracHits,
             'fracInteriorHits' : float(numInteriorHits)/float(numInteriorCells),
             'numConnectedInteriorCells' : numBothInboundOutbound,
             'fracConnectedInteriorCells' : float(numBothInboundOutbound)/float(numInteriorCells),
             'refMap' : refMap,
             'refMapId' : refMap.mapId}
        mapToData[refMap.mapId] = DictWrap(d)

    return mapToData

def dynamicProgammingTest2(queryMap, chunkDB, minDelta=1000, tol=0.1):
    """
    Compute the number of dynamic programming cells that are populated for a given query Map,
    after pruning.
    """

    # construct lower bounds and upper bounds
    frags = np.array(queryMap.frags[1:-1])
    numFrags = frags.shape[0]

    # Compute the upper/lower bounds for each fragment
    deltas = tol * frags
    minDeltas = minDelta * np.ones(numFrags)
    deltas = np.max([deltas, minDeltas], axis=0)
    lbs = frags - deltas
    ubs = frags + deltas

    fragHits = [chunkDB.getChunks(lb, ub) for lb,ub in zip(lbs, ubs)]

    # The coordinates (fragIndex, chromosome name, chromosome chunk end index) give the coordinates
    # of a single cell in the dynamic programming table.
    # For each chromosome, compute the total number of cells, and number of cells that have a match.
    coordGen = ((fInd, rc) for fInd, refHits in enumerate(fragHits) for rc in refHits)
    DPMap = lambda: defaultdict(list)
    refMapToHits = defaultdict(DPMap)
    for fInd, rc in coordGen:
        coord = (fInd, rc.eInd)
        destCoord = (fInd-1, rc.bInd)
        refMapToHits[rc.map][coord].append(destCoord)

    # Prune coord map
    def pruneDPMap(dpMap):
        done = False
        iterCount = 0
        while not done:
            madeChange = False
            DEBUG('pruneDpMap: iter=%i count=%i'%(iterCount, len(dpMap)))
            for coord, links in dpMap.items():
                # Only keep links that point to a cell that is in play or that
                # is beginning of query or reference
                # l[0]==-1 signifies beginning of query
                validLinks = [l for l in links if (l in dpMap) or (l[0]==-1)]
                madeChange = madeChange or (len(links) != len(validLinks))
                if validLinks:
                    dpMap[coord] = validLinks
                else:
                    del dpMap[coord]

            # Only keep cells that are pointed to.
            marked = set()
            for coord, links in dpMap.iteritems():
                marked.update(links)
            for coord in dpMap.keys():
                queryInd, refInd = coord
                # If coord in not pointed to and not in last row, remove it.
                if (coord not in marked) and (queryInd != numFrags-1):
                    del dpMap[coord]
                    madeChange = True
            done = not madeChange     
            iterCount += 1

    for map, dpMap in refMapToHits.iteritems():
        pruneDPMap(dpMap)

    mapToData = {}
    for map, hits in refMapToHits.iteritems():
        numMapFrags = map.numFrags
        numQueryFrags = numFrags
        numCells = numMapFrags*numQueryFrags
        numHits = len(hits)
        fracHits = float(numHits)/float(numCells)
        d = {'numMapFrags' : numMapFrags,
             'numQueryFrags' : numQueryFrags,
             'numCells' : numCells,
             'numHits' : numHits,
             'fracHits' : fracHits,
             'map' : map,
             'mapId' : map.mapId}
        mapToData[map] = DictWrap(d)

    return mapToData
