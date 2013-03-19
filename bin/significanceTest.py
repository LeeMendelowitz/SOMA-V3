##################################
# Author: Lee Mendelowitz (lmendelo@umiacs.umd.edu)
# Date: 1/22/2013
# File: significanceTest.py
#
# Description:
# Assign pvals to matches in the matchList by aligning random queries to the optical map
# and building a null distribution of random chunk scores.
# Random contigs are constructed by sampling with replacement from the distribution of contig restriction fragment lengths.
# Pvals for an alignment are computed by simulating scores for random alignments with the same number of alignment chunks 
# as the given alignment. The scores for a randomal alignment are simulated by sampling from the null distribution of random
# chunk scores.
#
#

import os
import numpy as np
import sys
import SOMAMap
import scipy
import subprocess
import parseSomaMatch

NUM_THREADS = 1
##############################################################################
# Class to store fragments
class FragDatabase(object):

    def __init__(self, *args, **kwargs):
        self.interiorFrags = np.array([])
        self.allFrags = np.array([])
        if 'mapFile' in kwargs:
            sys.stdout.write('Adding frags from map %s...'%kwargs['mapFile'])
            self.addFragsFromMap(kwargs['mapFile'])
            sys.stdout.write('Added %i fragments\n'%len(self.allFrags))

    def addInteriorFrags(self, frags):
        self.interiorFrags = np.array(list(self.interiorFrags) + list(frags))

    def addFrags(self, frags):
        self.allFrags = np.array(list(self.allFrags) + list(frags))

    # Sample n restriction fragments (with replacement)
    def sampleFrags(self, n):
        ind = np.random.randint(0, self.allFrags.shape[0], n)
        return self.allFrags[ind]
        
    # Sample n interior restriction fragments (with replacement)
    def sampleInteriorFrags(self, n):
        ind = np.random.randint(0, self.interiorFrags.shape[0], n)
        return self.interiorFrags[ind]

    def addFragsFromMap(self, mapFileName):
        mapDict = SOMAMap.readMaps(open(mapFileName))
        for mapId, mapObj in mapDict.iteritems():
            frags = mapObj.frags
            if not frags:
                continue
            self.addFrags(frags)
            if len(frags) > 2:
                self.addInteriorFrags(frags[1:-1])


################################################################################
# Class to store null distribution of contig chunk scores
class ChunkScoreNullDistribution(object):

    def __init__(self, chunkScores):
        self.chunkScores = np.sort(chunkScores)
        self.bootstrapScores = {}
        self.numSamples = 10000

    def computeScoreDistribution(self, numChunks):
        if numChunks in self.bootstrapScores:
            return self.bootstrapScores[numChunks]
        N = self.chunkScores.shape[0]
        scores = np.array( [np.sum(self.chunkScores[np.random.randint(0, N, numChunks)]) for i in xrange(self.numSamples)] )
        self.bootstrapScores[numChunks] = np.sort(scores)
        return scores

    # assign pval by computing how many random "alignments" built by sampling
    # N chunkScores from the empirical distribution of chunkScores score better than
    # this matchResult.
    def calcPval(self, mr):
        def assignPval(mr):
            score = mr.score
            numChunks = len(mr.alignment)
            scores = self.computeScoreDistribution(numChunks)
            numBetter = np.sum(scores >= score)
            pval = float(numBetter)/scores.shape[0]
            mr.pval = pval
            return pval

        if isinstance(mr, list):
            ml = mr
            for mr in ml:
                assignPval(mr)
        else:
            return assignPval(mr)

    def computeQuantile(self, numChunks, q):
        scores = self.computeScoreDistribution(numChunks)
        equantile = scipy.stats.mstats.mquantiles(scores, prob=[q])[0]
        return equantile

###########################################################
def generateRandomMaps(mapFileIn, numRandomMaps, fragsPerMap, mapFileOut):
    fout = open(mapFileOut, 'w')
    fragDB = FragDatabase(mapFile=mapFileIn)

    for i in range(numRandomMaps):
        frags = list(fragDB.sampleInteriorFrags(fragsPerMap))
        # Bound by fragments of length 0, indicating that this random contig starts and ends with a restriction site
        #frags = [0] + frags  + [0]
        mapId = 'map%i_frags%i'%(i, fragsPerMap)
        map = SOMAMap.SOMAMap(mapId=mapId, frags=frags)
        map.write(fout)
    fout.close()

###########################################################
# Wrapper around SOMA to align a random query map to an optical map.
# The SOMA settings are set to only take the best match for a random
# query map. All restriction fragments are treated as interior restriction
# fragments (bounded by two restriction sites)
def alignRandomMaps(randomQueryMapFile, opticalMapFile, outputPfx, numThreads = NUM_THREADS):
    cmd = ['./match',
           '--maxMatchesPerContig=1',
           '--siteCostContig=3',
           '--siteCostOptical=5',
           '--sdMin=4',
           '--sdMax=4',
           '--noBoundaries',
           '--maxMissRateContig=1.00',
           '--minSiteHits=1',
           '--numThreads=%i'%numThreads,
           '--output=%s'%outputPfx]
    cmd.append(randomQueryMapFile)
    cmd.append(opticalMapFile)
    sys.stderr.write('CMD: %s\n'%(' '.join(cmd)))
    subprocess.call(cmd)

#############################################################
# Parse SOMA Alignments.
# Return a dictionary which stores:
#   matchDict: a dict with keys=queryMapId, value = match result
#   chunkScores: a list of the scores if alignment chunks
#   scores: a list of the total score for each alignment
# Note: This function expects there to be only one match per query
def parseRandomAlignments(xmlFile):
    numFrags = int(xmlFile.split('.')[1])
    matchResults = parseSomaMatch.parseMatchFileXML(xmlFile)
    ml = matchResults
    matchDict = dict( (mr.contigId, mr) for mr in matchResults)
    resDict = {}
    resDict['matchDict'] = matchDict
    resDict['chunkScores'] = [chunkScore for mr in ml for chunkScore in mr.getChunkScores()]
    resDict['scores'] = [mr.score for mr in ml]
    return resDict

#############################################################
# Build the empirical distribution of chunk scores for random
# alignments by generating random query maps and collecting the
# scores for each alignment chunk.
# Return a ChunkScoreNullDistribution instance
def buildChunkScoreNullDistribution(contigMapFile, opticalMapFile, numRandomContigs, numFragsPerContig):
    randomContigMapFile = 'randomContigs.%i.silico'%numFragsPerContig
    alnPfx = '%s.%s'%(randomContigMapFile, opticalMapFile)
    xmlFile = '%s_AllMatches.xml'%alnPfx

    if os.path.exists(xmlFile):
        sys.stderr.write('XML File for random alignments already exists: %s. Using these alignments.\n'%xmlFile)
    else:
        generateRandomMaps(contigMapFile, numRandomContigs, numFragsPerContig, randomContigMapFile)
        alignRandomMaps(randomContigMapFile, opticalMapFile, alnPfx)
    randomAlnData = parseRandomAlignments(xmlFile)
    return ChunkScoreNullDistribution(randomAlnData['chunkScores'])

#############################################################
# Assign pvals to matches in the matchList by aligning random queries to the optical map
# and building a null distribution of random chunk scores.
# Random contigs are constructed by sampling with replacement from the distribution of contig restriction fragment lengths.
# Pvals for an alignment are computed by simulating scores for random alignments with the same number of alignment chunks 
# as the given alignment. The scores for a randomal alignment are simulated by sampling from the null distribution of random
# chunk scores.
def runSignificanceTest(matchList, contigMapFile, opticalMapFile, numRandomContigs=1000, numFragsPerContig=10, numThreads = 1):
    global NUM_THREADS
    NUM_THREADS = numThreads
    nullDist = buildChunkScoreNullDistribution(contigMapFile, opticalMapFile, numRandomContigs, numFragsPerContig)
    nullDist.calcPval(matchList)
