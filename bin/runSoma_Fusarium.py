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
#from ruffus import *
from collections import Counter

# For parsing output
import parseSomaMatch
import summarizeContigStatus
import createScaffolds
import make_opt
import make_silico
import SOMAMap
import significanceTest
import alignRandomMaps
import adjustOpticalMaps

CWD = os.getcwd()
SOMA_BIN_DIR = CWD

# SOMA Options
numThreads = 4

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


###################################################
# Convert an optical map from the Schwartz lab format
# to the SOMA format
# opticalMapFile: optical map file in the Schwartz lab format
def convertOpticalMaps(opticalMapFile, outputPfx):
    opMapFileOut = '%s.opt'%outputPfx

    msg = '\n'+'*'*50 + \
          '\nReading Optical Map File %s\n'%opticalMapFile + \
          '*'*50 + '\n'
    sys.stderr.write(msg)

    opMapList = make_opt.readMapDataSchwartz(opticalMapFile)
    enzymeSet = set(om.enzyme for om in opMapList)
    if len(enzymeSet) > 1:
        raise RuntimeError('Different enzymes used in the input optical map set!')
    enzyme = opMapList[0].enzyme

    msg = '\n'+'*'*50 +\
          '\nConverting Optical Map to SOMA Format\n' +\
          '*'*50 + '\n'
    sys.stderr.write(msg)

    # Optical maps for chromosomes 
    # Remove all white space from restriction map names
    for opMap in opMapList:
        opMap.mapId = ''.join(opMap.mapId.split())
    SOMAMap.writeMaps(opMapList, opMapFileOut)
    result = { 'enzyme' : enzyme,
               'opMapList' : opMapList,
               'opticalMapFile' : opMapFileOut}
    return result


#############################
# opticalMapFile: optical map file in the Schwartz lab format
# contigFile: fasta with contigs to be aligned to the optical map
def createMaps(opticalMapFileIn, contigFile, outputPfx):
    contigMapFile = '%s.silico'%outputPfx
    res = convertOpticalMaps(opticalMapFileIn, outputPfx)
    opticalMapFile = res['opticalMapFile']
    enzyme = res['enzyme']
    if not os.path.exists(contigMapFile):
        msg = '\n'+'*'*50+'\n' + \
              'Creating in-silico optical map using contig file %s and enzyme %s\n'%(contigFile, enzyme)+\
              '*'*50 + '\n'
        sys.stderr.write(msg)
        make_silico.makeInsilicoMain(contigFile, enzyme, contigMapFile)
    else:
        msg = 'Warning: Contig map file %s already exists. Using this map file.\n'
        sys.stderr.write(msg)
    result = { 'contigMapFile': contigMapFile,
               'opticalMapFile': opticalMapFile }
    return result

########################################################
# Make alignments of the maps in the contigMap file to the maps in the opticalMap file
# opticalMap: file in SOMA Map format
# contigMap: map file in SOMA Map format
def makeAlignments(opticalMap, contigMap, outputPfx):

    ###########################################################################
    xmlFile = '%s_SigMatches.xml'%outputPfx
    result = { 'xmlFile': xmlFile}

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

    return result

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
    significanceTest.runSignificanceTest(ml, contigMapFile, opticalMapFile, numThreads=numThreads)

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
    res = createMaps(opticalMapFile, contigFile, pfx)
    opMapSoma = res['opticalMapFile']
    cMapSoma = res['contigMapFile']
    makeAlignments(opMapSoma, cMapSoma, pfx)

    xmlFile = '%s_AllMatches.xml'%pfx

    postProcess(xmlFile, opMapSoma, cMapSoma, pfx)

def runSomaFusarium_2_23():
    wd = '/cbcb/project-scratch/lmendelo/Fusarium/broadAssembly/hybridAssembly_Oct/output7/selfNucmerAlignment/newScaffs/somaAlignment'
    os.chdir(wd)
    opticalMapFile =  'Fusarium_oxysporiumOpticalMaps061306'
    contigSeqFile = 'newScaffolds.fasta'
    pfx = 'aln1'

    # Convert into SOMA formats
    res = createMaps(opticalMapFile, contigSeqFile, pfx)
    opMapSoma = res['opticalMapFile']
    cMapSoma = res['contigMapFile']

    # Adjust the optical map
    opMapSomaAdjusted = '%s.adjusted.opt'%pfx
    adjustOpticalMap(cMapSoma, opMapSoma, opMapSomaAdjusted)

    # Make alignments against the adjusted map
    res = makeAlignments(opMapSomaAdjusted, cMapSoma, pfx)
    xmlFile = res['xmlFile']
    postProcess(xmlFile, opMapSoma, cMapSoma, pfx)

def runSomaFusarium_2_26():
    wd = '/cbcb/project-scratch/lmendelo/Fusarium/broadAssembly/hybridAssembly_Oct/output7/selfNucmerAlignment/newScaffs/reScaffold/somaAlignment'
    os.chdir(wd)
    opticalMapFile =  'Fusarium_oxysporiumOpticalMaps061306'
    contigSeqFile = 'reScaffolds.fasta'
    pfx = 'aln_reScaffolds'

    # Convert into SOMA formats
    res = createMaps(opticalMapFile, contigSeqFile, pfx)
    opMapSoma = res['opticalMapFile']
    cMapSoma = res['contigMapFile']

    # Adjust the optical map
    opMapSomaAdjusted = '%s.adjusted.opt'%pfx
    adjustOpticalMap(cMapSoma, opMapSoma, opMapSomaAdjusted)

    # Make alignments against the adjusted map
    res = makeAlignments(opMapSomaAdjusted, cMapSoma, pfx)
    xmlFile = res['xmlFile']
    postProcess(xmlFile, opMapSoma, cMapSoma, pfx)

def runSomaFusarium_3_10():
    wd = '/cbcb/project-scratch/lmendelo/Fusarium/broadAssembly/hybridAssembly_Oct/output7/contigCleanup/bambusMergedContigScaffolds/opMapAlign'
    os.chdir(wd)
    opticalMapFile =  'Fusarium_oxysporiumOpticalMaps061306'
    contigSeqFile = 'bambusScaff.scaffolds.linear.formatted.fasta'
    pfx = 'aln_bambusScaffolds'

    # Convert into SOMA formats
    res = createMaps(opticalMapFile, contigSeqFile, pfx)
    opMapSoma = res['opticalMapFile']
    cMapSoma = res['contigMapFile']

    # Adjust the optical map
    opMapSomaAdjusted = '%s.adjusted.opt'%pfx
    adjustOpticalMap(cMapSoma, opMapSoma, opMapSomaAdjusted)

    # Make alignments against the adjusted map
    res = makeAlignments(opMapSomaAdjusted, cMapSoma, pfx)
    xmlFile = res['xmlFile']
    postProcess(xmlFile, opMapSoma, cMapSoma, pfx)

def runAdjustedMaps():
    opMapSoma = 'aln0_global.adjust.opt' 
    cMapSoma = 'aln0_global.silico'
    pfx = 'aln0_global_adjusted'
    xmlFile = '%s_AllMatches.xml'%pfx
    makeAlignments(opMapSoma, cMapSoma, pfx)
    postProcess(xmlFile, opMapSoma, cMapSoma, pfx)

###############################################################
# Perform an initial alignment of the insilico maps
# Use high quality alignments 
# Use the methods in adjustOpticalMaps module to create a new optical map
def adjustOpticalMap(contigMapFile, opticalMapFileOrig, opticalMapFileNew):

    # Perform an initial alignment of of the contigMapFile against the opticalMapFile
    outputPfx = '%s.%s.alignForAdjustment'%(contigMapFile,opticalMapFileOrig)
    res = makeAlignments(opticalMapFileOrig, contigMapFile, outputPfx)
    ml = parseSomaMatch.parseMatchFileXML(res['xmlFile'])

    sys.stderr.write('adjustOpticalMap: Parsed %i matches from file %s\n'%(len(ml), res['xmlFile']))

    # Filter the match list:
    minHits = 10
    maxMissRate = 0.10
    def matchOK(mr):
        hitsOK = mr.contigHits >= minHits
        missRateOK = max(mr.contigMissRate, mr.opticalMissRate) <= maxMissRate
        return hitsOK and missRateOK
    goodMatches = [mr for mr in ml if matchOK(mr)]
    sys.stderr.write('adjustOpticalMap: Filtered to %i matches based on quality'%(len(goodMatches)))

    # Count the number of good alignments per contig. Only select the unique alignments
    contigAlignmentCounts = Counter(mr.contigId for mr in goodMatches)
    uniqueAlignments = [mr for mr in goodMatches if contigAlignmentCounts[mr.contigId]==1]
    sys.stderr.write('adjustOpticalMap: Filtered to %i matches based on uniqueness'%(len(uniqueAlignments)))

    # Create a matched chunk file
    matchedChunkFile = '%s.matchedChunks'%outputPfx
    adjustOpticalMaps.makeMatchedChunkFile(uniqueAlignments, matchedChunkFile)
    adjustOpticalMaps.run(opticalMapFileOrig, matchedChunkFile, opticalMapFileNew)

# align optical map to itself
def selfAlignment():
    outputPfx = 'opMapSelfAlignment'
    opticalMapFile = 'Fusarium_oxysporiumOpticalMaps061306'
    opticalMapFileSOMA = '%s.opt'%outputPfx
    xmlFile = '%s_AllMatches.xml'%outputPfx

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
    SOMAMap.writeMaps(opMapList, opticalMapFileSOMA)

    # Run the self alignment
    global local
    local = True
    makeAlignments(opticalMapFileSOMA, opticalMapFileSOMA, outputPfx)
    postProcess(xmlFile, opticalMapFileSOMA, opticalMapFileSOMA, outputPfx)
    
def main():
    runSomaFusarium_3_10()
    return
    #opticalMap, contigFile, outputPfx = sys.argv[1:4]
    opticalMap =  'Fusarium_oxysporiumOpticalMaps061306'
    #contigFile =  'fusarium.oxysporum4287.v1p3.contigs.fa'
    contigFile = 'fusarium.oxysporum4287.v1p3.scaffolds.fa'
    outputPfx = 'aln0_global'
    runSoma(opticalMap, contigFile, outputPfx)

if __name__ == '__main__':
   main()
