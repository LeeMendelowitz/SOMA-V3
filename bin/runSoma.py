#! /usr/bin/env python

# SOMA Match Pipeline
# Inputs: Optical Map File (Schwartz Lab Format) and Sequence Fasta File
# This will convert the optical map into the Soma Format
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

CWD = os.getcwd()
SOMA_BIN_DIR = CWD

# SOMA Options
numThreads = 1
numPermutationTrials = 0

# Parameters
sdMin=4
sdMax=4
maxMissRateContig = 0.50
pThreshold = 0.05
siteCostContig = 3
siteCostOptical = 5

def getBaseName(inFile):
    (dirName, fileName) = os.path.split(inFile)
    (fileBaseName, fileExtension)=os.path.splitext(fileName)
    return fileBaseName

def runSoma(opticalMapFile, contigFile, pfx):

    silicoFileOut = '%s.silico'%pfx
    opticalMapFileOut = '%s.opt'%getBaseName(opticalMapFile)

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
    SOMAMap.writeMaps(opMapList, opticalMapFileOut)
    ###########################################################################
    print '\n'+'*'*50
    print 'Creating in-silico optical map using contig file %s and enzyme %s'%(contigFile, enzyme)
    print '*'*50 + '\n'

    make_silico.makeInsilicoMain(contigFile, enzyme, silicoFileOut)

    ###########################################################################
    print '\n'+'*'*50
    print 'Matching contigs to the optical maps'
    print '*'*50 + '\n'

    # Add commands to map to individual optical maps
    outputPfx = pfx
    cmd = [os.path.join(SOMA_BIN_DIR, 'match'),
    '--output', outputPfx,
    '--sdMin', str(sdMin),
    '--sdMax', str(sdMax),
    '--pvalue', str(pThreshold),
    '--maxMissRateContig', str(maxMissRateContig),
    '--numThreads', str(numThreads),
    '--numPermutationTrials',str(numPermutationTrials),
    '--siteCostContig', str(siteCostContig),
    '--siteCostOptical', str(siteCostOptical)]
    cmd.append(silicoFileOut)
    cmd.append(opticalMapFileOut)
    print 'CMD: ',' '.join(cmd)
    retCode = subprocess.call(cmd)

    if retCode != 0:
        raise RuntimeError('SOMA MATCH returned with code %i'%retCode)

    ##################################################################################
    # Parse results
    print '\n'+'*'*50
    print 'Parsing SOMA OUTPUT...'
    print '*'*50 + '\n'

    bn = outputPfx + '_SigMatches'
    xmlFile = bn + '.xml'
    pickleFileUnique = bn + '.refinedUniqueMatchList.pickle'
    pickleFileAll = bn + '.refinedMatchList.pickle'

    # Parse Match Results. Write Pickle Files
    parseSomaMatch.parse(xmlFile)

    # Load Match Lists
    mlUnique = cPickle.load(open(pickleFileUnique))
    mlAll = cPickle.load(open(pickleFileAll))

    # Summarize alignment status for contigs in the silicoFile
    summarizeContigStatus.summarizeContigStatus(bn, silicoFileOut)

    #  Print all of the alignments to a textFile
    fout = open(bn + '.UniqueAlignments.txt', 'w')
    parseSomaMatch.printAlignments(mlUnique, fout)
    fout.close()
    fout = open(bn + '.AllAlignments.txt', 'w')
    parseSomaMatch.printAlignments(mlAll, fout)
    fout.close()

    # Create scaffolds
    print '\n'+'*'*50
    print 'Creating Scaffolds...'
    print '*'*50 + '\n'
    createScaffolds.createScaffolds(pickleFileAll, opMapFiles, '%s.scaffold_allAlignments_withOverlaps.txt'%bn, allowOverlaps=True)
    createScaffolds.createScaffolds(pickleFileAll, opMapFiles, '%s.scaffold_allAlignments_noOverlaps.txt'%bn, allowOverlaps=False)
    createScaffolds.createScaffolds(pickleFileUnique, opMapFiles, '%s.scaffold_UniqueAlignments_withOverlaps.txt'%bn, allowOverlaps=True)
    createScaffolds.createScaffolds(pickleFileUnique, opMapFiles, '%s.scaffold_UniqueAlignments_noOverlaps.txt'%bn, allowOverlaps=False)


def main():
    opticalMap, contigFile, outputPfx = sys.argv[1:4]
    runSoma(opticalMap, contigFile, outputPfx)

if __name__ == '__main__':
   main()
