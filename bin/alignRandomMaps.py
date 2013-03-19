# Generate random contig maps by sampling from the
# distribution of interior contig fragments
import subprocess
import sys
import matplotlib.pyplot as plt
from significanceTest import FragDatabase
from significanceTest import ChunkScoreNullDistribution
import SOMAMap
import MatchResult
import parseSomaMatch
import numpy as np
import os

NUM_THREADS=4

###########################################################
def generateRandomMaps(mapFileIn, numMaps, fragsPerMap, mapFileOut):
    fout = open(mapFileOut, 'w')
    fragDB = FragDatabase(mapFile=mapFileIn)

    for i in range(numMaps):
        frags = list(fragDB.sampleInteriorFrags(fragsPerMap))
        # Bound by fragments of length 0, indicating that this random contig starts and ends with a restriction site
        #frags = [0] + frags  + [0]
        mapId = 'map%i_frags%i'%(i, fragsPerMap)
        map = SOMAMap.SOMAMap(mapId=mapId, frags=frags)
        map.write(fout)
    fout.close()

###########################################################
# Run SOMA to align the random contig maps. Only accept 
# the best match per random map, construct the match result
def alignRandomMaps(randomContigFile, opticalMapFile, outputPfx):
    cmd = ['./match',
           '--maxMatchesPerContig=1',
           '--siteCostContig=3',
           '--siteCostOptical=5',
           '--sdMin=4',
           '--sdMax=4',
           '--noBoundaries',
           '--maxMissRateContig=1.00',
           '--minSiteHits=1',
           '--numThreads=%i'%NUM_THREADS,
           '--output=%s'%outputPfx]
    cmd.append(randomContigFile)
    cmd.append(opticalMapFile)
    sys.stderr.write('CMD: %s\n'%(' '.join(cmd)))
    subprocess.call(cmd)

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
# alignments by generating random contigs and collected the
# scores of MatchedChunks
def buildChunkScoreNullDistribution(contigMapFile, opticalMapFile, numRandomContigs, numFragsPerContig):
    randomContigMapFile = 'randomContigs.%i.silico'%numFragsPerContig
    alnPfx = '%s.%s'%(randomContigMapFile, opticalMapFile)
    xmlFile = '%s_AllMatches.xml'%alnPfx

    if os.path.exists(xmlFile):
        sys.stderr.write('XML File for random alignments already exists: %s. Using these alignments.\n'%xmlFile)
    else:
        generateRandomMaps(contigMapFile, numRandomContigs, numFragsPerContig, randomContigMapFile)
        alignRandomMaps(randomContigMapFile, opticalMapFile, alnPfx)

    ml = parseSomaMatch.parseMatchFileXML(xmlFile)
    chunkScores = [chunkScore for mr in ml
                              for chunkScore in mr.getChunkScores()]
    chunkScores = np.sort(chunkScores)
    return ChunkScoreNullDistribution(chunkScores)


def run1():
    contigMapFile = 'aln0_global.silico'
    numMaps = 1000
    fragsPerMap = 5
    randomContigMapFile = 'randomContigs.%i.silico'%fragsPerMap
    alnPfx = 'alignRandom'
    opticalMapFile = 'aln0_global.opt'
    generateRandomMaps(contigMapFile, numMaps, fragsPerMap, randomContigMapFile)
    matchDict = alignRandomMaps(randomContigMapFile, opticalMapFile, alnPfx)
    return matchDict

def run2():
    contigMapFile = 'aln0_global.silico'
    opticalMapFile = 'aln0_global.opt'
    numMaps = 1000
    fragsPerMap = [5,10,15,20,25,30,40,50]
    fragsPerMap = [5,15]
    results = []
    for numFrags in fragsPerMap:
        resDict = {}
        randomContigMapFile = 'randomContigs.%i.silico'%numFrags
        alnPfx = '%s-%s'%(randomContigMapFile, opticalMapFile)
        generateRandomMaps(contigMapFile, numMaps, numFrags, randomContigMapFile)
        alignRandomMaps(randomContigMapFile, opticalMapFile, alnPfx)
        xmlFile = '%s_AllMatches.xml'%alnPfx
        resDict = parseRandomAlignments(xmlFile)
        resDict['numFrags'] = numFrags
        resDict['contigMapFile'] = randomContigMapFile
        results.append(resDict)
    return results

#####################################################################################
# Assign a pval to matchResults using a null distribution constructed from random alignments
def run3():
    contigMapFile = 'aln0_global.silico'
    opticalMapFile = 'aln0_global.opt'
    randomAlignmentXml = 'randomContigs.20.silico-aln0_global.opt_AllMatches.xml'
    #randomAlignmentXml = 'randomContigs.5.silico-aln0_global.opt_AllMatches.xml'
    contigAlignmentXml = 'aln0_global_AllMatches.xml'

    # Build null distriubtion
    randomAlns = parseRandomAlignments(randomAlignmentXml)
    nullDist = ChunkScoreNullDistribution(randomAlns['chunkScores'])

    # Read contig alignments
    ml = parseSomaMatch.parseMatchFileXML(contigAlignmentXml)
    for mr in ml:
        nullDist.calcPval(mr)
    return ml


# Plot a pval of 0.05 cutoff using the chunk score null distribution
# sampled from different random contigs
def plotCutoffs():
    contigMapFile = 'aln0_global.silico'
    opticalMapFile = 'aln0_global.opt'
    randomAlignmentXmlTemplate = 'randomContigs.%i.silico-aln0_global.opt_AllMatches.xml'

    # Build null distriubtion
    ncList = [5, 10, 15, 20, 25]
    legendLabels = []
    plt.figure()
    for numContigChunks in ncList:
        randomAlignmentXml = randomAlignmentXmlTemplate%numContigChunks
        randomAlns = parseRandomAlignments(randomAlignmentXml)
        nullDist = ChunkScoreNullDistribution(randomAlns['chunkScores'])
        numChunks = np.arange(5, 301,5)
        pval = 0.05
        cutoff = [nullDist.computeQuantile(nc, 1-pval) for nc in numChunks]
        #plt.figure()
        label = '%i_chunks'%numContigChunks
        plt.plot(numChunks, cutoff, label=label)
        legendLabels.append(label)
    plt.title('Cutoff scores, pval = %f'%pval)
    plt.xlabel('Num Alignment Chunks')
    plt.legend(tuple(legendLabels))
    plt.show()



def plotScoreBoxPlot(resList):
    plt.figure()
    plt.title('Score Distributions for random contigs')
    scoreVec = [res['scores'] for res in resList]
    numFrags = [res['numFrags'] for res in resList]
    xlabels = [str(n) for n in numFrags]
    N = len(scoreVec)
    plt.boxplot(scoreVec)
    plt.xticks(range(1,N+1), xlabels)
    plt.xlabel('Num Fragments per contig')

    plt.figure()
    plt.title('Chunk Score Distributions for random contigs')
    scoreVec = [res['chunkScores'] for res in resList]
    numFrags = [res['numFrags'] for res in resList]
    xlabels = [str(n) for n in numFrags]
    N = len(scoreVec)
    plt.boxplot(scoreVec)
    plt.xticks(range(1,N+1), xlabels)
    plt.xlabel('Num Fragments per contig')




if __name__ == '__main__':
    run1()
