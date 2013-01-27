#! /usr/bin/env python

# SOMA Match Pipeline
# Inputs: Optical Map File (Schwartz Lab Format) and Sequence Fasta File
# This will convert the optical map into the SOMA Format
# and align contigs/scaffolds from the sequence file to the optical map.

import subprocess
import os
import glob
import sys
import cPickle
from ruffus import *

# For parsing output
import parseSomaMatch
import summarizeContigStatus
import createScaffolds
import make_opt
import make_silico
import SOMAMap
import alignRandomMaps

CWD = os.getcwd()
SOMA_BIN_DIR = CWD

# SOMA Options
numThreads = 4
numPermutationTrials = 0

# Parameters
sdMin=4
sdMax=4
maxMissRateContig = 0.50
maxMatchesPerContig = 5
minSiteHits  = 2
pThreshold = 0.00
siteCostContig = 3
siteCostOptical = 5
local = False

def getBaseName(inFile):
    (dirName, fileName) = os.path.split(inFile)
    (fileBaseName, fileExtension)=os.path.splitext(fileName)
    return fileBaseName

###########################
# Assign pvals to each match result in the matchlist
def assignPvals(ml, contigMapFile, opticalMapFile):
    fragsPerRandomContig = 10
    numRandomContigs = 1000
    nullDist = alignRandomMaps.buildChunkScoreNullDistribution(contigMapFile, opticalMapFile, numRandomContigs, fragsPerRandomContig)
    nullDist.calcPval(ml)

#############################
# opticalMapFile: optical map file in the Schwartz lab format
# contigFile: fasta with contigs to be aligned to the optical map
def createMaps(opticalMapFile, contigFile, outputPfx):
    silicoFileOut = '%s.silico'%outputPfx
    opMapFileOut = '%s.opt'%outputPfx

    ###########################################################################
    print '\n'+'*'*50
    print 'Reading Optical Map File %s'%opticalMapFile
    print '*'*50 + '\n'

    opMapList = make_opt.readMapDataSchwartz(opticalMapFile)

    # Assert that all of the optical maps have the same enzyme list.
    # otherwise this pipeline will not work
    assert len(set([om.enzyme for om in opMapList])) == 1
    enzyme = opMapList[0].enzyme

    ##########################################################################
    print '\n'+'*'*50
    print 'Converting Optical Map to SOMA Format'
    print '*'*50 + '\n'

    # Optical maps for chromosomes 
    # Remove all white space from restriction map names
    for opMap in opMapList:
        oldName = opMap.mapId
        opMap.mapId = ''.join(oldName.split())
    SOMAMap.writeMaps(opMapList, opMapFileOut)

    ###########################################################################
    print '\n'+'*'*50
    print 'Creating in-silico optical map using contig file %s and enzyme %s'%(contigFile, enzyme)
    print '*'*50 + '\n'
    if not os.path.exists(silicoFileOut):
        make_silico.makeInsilicoMain(contigFile, enzyme, silicoFileOut)


########################################################
# Make alignments of the maps in the contigMap file to the maps in the opticalMap file
# opticalMap: file in SOMA Map format
# contigMap: map file in SOMA Map format
def makeAlignments(opticalMap, contigMap, outputPfx):

    silicoFileOut = '%s.silico'%outputPfx
    opMapFileOut = '%s.opt'%outputPfx


    ###########################################################################
    xmlFile = '%s_SigMatches.xml'%outputPfx

    print '\n'+'*'*50
    print 'Matching contigs to the optical maps'
    print '*'*50 + '\n'

    # Add commands to map to individual optical maps
    cmd = [os.path.join(SOMA_BIN_DIR, 'match'),
    '--output', outputPfx,
    '--sdMin', str(sdMin),
    '--sdMax', str(sdMax),
    '--pvalue', str(pThreshold),
    '--maxMissRateContig', str(maxMissRateContig),
    '--maxMatchesPerContig=%i'%maxMatchesPerContig,
    '--minSiteHits=%i'%minSiteHits,
    '--numThreads', str(numThreads),
    '--numPermutationTrials',str(numPermutationTrials),
    '--siteCostContig', str(siteCostContig),
    '--siteCostOptical', str(siteCostOptical)]

    if local:
        cmd.append('--local')

    cmd.append(contigMap)
    cmd.append(opticalMap)

    print 'CMD: ',' '.join(cmd)
    if not os.path.exists(xmlFile):
        retCode = subprocess.call(cmd)

        if retCode != 0:
            raise RuntimeError('SOMA MATCH returned with code %i'%retCode)
    else:
        sys.stderr.write('Output file %s already exists. Not running SOMA.\n'%xmlFile)

############################################################################
def postProcess(xmlFile, opticalMapFile, contigMapFile, outputPfx):

    # Parse results
    print '\n'+'*'*50
    print 'Parsing SOMA OUTPUT...'
    print '*'*50 + '\n'

    pickleFileAll = '%s.matchList.all.pickle'%outputPfx
    pickleFileSig = '%s.matchList.sig.pickle'%outputPfx
    pickleFileSigUnique = '%s.matchList.sig.unique.pickle'%outputPfx

    # Parse Match Results. Write Pickle Files
    ml = parseSomaMatch.parseMatchFileXML(xmlFile)
    assignPvals(ml, contigMapFile, opticalMapFile)

    # Select significant results
    pvalCutoff = 0.05
    sigMatches = [mr for mr in ml if mr.pval <= pvalCutoff]
    sigMatchDict = parseSomaMatch.collectMatchResultsByContig(sigMatches)
    sigUniqueMatches = [matches[0] for contigId, matches in sigMatchDict.iteritems() if len(matches)==1]
    sigUniqueMatchDict = parseSomaMatch.collectMatchResultsByContig(sigUniqueMatches)
    sys.stdout.write('Found %i significant matches (%i bp)\n'%(len(sigMatches), sum(mr.cAlignedBases for mr in sigMatches)))
    sys.stdout.write('Found %i unique significant matches (%i bp)\n'%(len(sigUniqueMatches), sum(mr.cAlignedBases for mr in sigUniqueMatches)))

    # Pickle the matchResults
    cPickle.dump(ml, open(pickleFileAll, 'w'))
    cPickle.dump(sigMatches, open(pickleFileSig, 'w'))
    cPickle.dump(sigUniqueMatches, open(pickleFileSigUnique, 'w'))

    infoFileOut = '%s.info'%outputPfx
    parseSomaMatch.writeInfoFile2(ml, infoFileOut)

    # Summarize alignment status for contigs in the silicoFile
    contigMapDict = SOMAMap.readMaps(contigMapFile)
    opMapDict= SOMAMap.readMaps(opticalMapFile)
    summarizeContigStatus.summarizeContigStatus(outputPfx, sigMatchDict, contigMapDict)

    #  Print all of the alignments to a textFile
    fout = open('%s.SigUniqueAlignments.txt'%outputPfx, 'w')
    parseSomaMatch.printAlignments(sigUniqueMatches, fout)
    fout.close()
    fout = open('%s.AllSigAlignments.txt'%outputPfx, 'w')
    parseSomaMatch.printAlignments(sigMatches, fout)
    fout.close()

    # Create scaffolds
    print '\n'+'*'*50
    print 'Creating Scaffolds...'
    print '*'*50 + '\n'

    createScaffolds.createScaffolds(sigMatchDict, opMapDict, '%s.scaffold_sigMatches_withOverlaps.txt'%outputPfx, allowOverlaps=True)
    createScaffolds.createScaffolds(sigMatchDict, opMapDict, '%s.scaffold_sigMatches_noOverlaps.txt'%outputPfx, allowOverlaps=False)
    createScaffolds.createScaffolds(sigUniqueMatchDict, opMapDict, '%s.scaffold_sigUniqueMatches_withOverlaps.txt'%outputPfx, allowOverlaps=True)
    createScaffolds.createScaffolds(sigUniqueMatchDict, opMapDict, '%s.scaffold_sigUniqueMatches_noOverlaps.txt'%outputPfx, allowOverlaps=False)

def runSoma(opticalMapFile, contigFile, pfx):
    createMaps(opticalMapFile, contigFile, pfx)

    opMapSoma = '%s.opt'%pfx
    cMapSoma = '%s.silico'%pfx
    xmlFile = '%s_AllMatches.xml'%pfx
    makeAlignments(opMapSoma, cMapSoma, pfx)

    postProcess(xmlFile, opMapSoma, cMapSoma, pfx)

def runAdjustedMaps():
    opMapSoma = 'aln0_global.adjust.opt' 
    cMapSoma = 'aln0_global.silico'
    pfx = 'aln0_global_adjusted'
    xmlFile = '%s_AllMatches.xml'%pfx
    makeAlignments(opMapSoma, cMapSoma, pfx)
    postProcess(xmlFile, opMapSoma, cMapSoma, pfx)


def main():
    #opticalMap, contigFile, outputPfx = sys.argv[1:4]
    opticalMap =  'Fusarium_oxysporiumOpticalMaps061306'
    #contigFile =  'fusarium.oxysporum4287.v1p3.contigs.fa'
    contigFile = 'fusarium.oxysporum4287.v1p3.scaffolds.fa'
    outputPfx = 'aln0_global'
    runSoma(opticalMap, contigFile, outputPfx)

if __name__ == '__main__':
   main()
