########################################################################
# Filename: parseSomaMatch.py
# Author: Lee Mendelowitz
# Date: 4/8/2012

# Script for parsing SOMA's .unique_match or .all_match output file
# Creates 
########################################################################
import numpy as np
import sys
import glob
import re
import os
import parallelProcess as pp
import cPickle
from operator import attrgetter
import matplotlib.pyplot as plt
import xml.etree.cElementTree as ElementTree
from MatchResult import MatchResult
import importSomaData

########################################################################
def getBaseName(inFile):
    (dirName, fileName) = os.path.split(inFile)
    (fileBaseName, fileExtension)=os.path.splitext(fileName)
    return fileBaseName

########################################################################
def parseMatchFileXMLQuality(filename):
    matchResultList = []
    matchCount = 0
    # Use the XML parser to parse the file. To save memory, delete
    # the xmlElement for the MatchResult once the MatchResult has been
    # instantiated
    for event,element in ElementTree.iterparse(filename):
        if element.tag == 'MatchResult':
            mr = MatchResult(xmlElement = element)
            element.clear() # Clear the xml Element
            matchCount += 1
            if isQualityMatch(mr):
                matchResultList.append(mr)
    print 'File: %s. Read %i matches, found %i quality matches'%(filename, matchCount, len(matchResultList))
    return matchResultList

########################################################################
def parseMatchFileXML(filename):
    matchResultList = []
    matchCount = 0
    # Use the XML parser to parse the file. To save memory, delete
    # the xmlElement for the MatchResult once the MatchResult has been
    # instantiated
    matchResultList = [MatchResult(xmlElement=element) for event,element in ElementTree.iterparse(filename) if element.tag == 'MatchResult']
    '''
    for event,element in ElementTree.iterparse(filename):
        if element.tag == 'MatchResult':
            mr = MatchResult(xmlElement = element, chromosome=chromName)
            element.clear() # Clear the xml Element
            matchCount += 1
            matchResultList.append(mr)
    '''
    matchCount = len(matchResultList)
    print 'File: %s. Read %i matches.'%(filename, matchCount)
    return matchResultList

LENGTH_RATIO = 0.90 # A MatchResults length ratio must be greater than or equal to this
MISS_RATE = 1.0 # This fraction of missed sites must be less than this
UNALIGNED_RATIO = 0.10 # A MatchResult fraction of contig bases that are unaligned must be less than this
PVAL = 0.05

########################################################################
# Return True if the MatchResult mr is a quality match
# False otherwise
def isQualityMatch(mr):
    if mr.numAlignedFrags<3:
        return False
    lengthCheck = mr.lengthRatio >= LENGTH_RATIO
    missRateCheck = mr.totalMissRate <= MISS_RATE
    pvalCheck = mr.pval <= PVAL
    unAlignedCheck = mr.contigUnalignedBaseRatio <= UNALIGNED_RATIO
    return (lengthCheck and missRateCheck and unAlignedCheck and pvalCheck)

########################################################################
# Stricter requirements on oneToOneCheck
def isQualityMatchA(mr):
    if mr.numAlignedFrags<3:
        return False
    lengthCheck = mr.lengthRatio >= LENGTH_RATIO
    missRateCheck = mr.totalMissRate <= MISS_RATE
    unAlignedCheck = mr.contigUnalignedBaseRatio <= UNALIGNED_RATIO
    oneToOneCheck = mr.basesInOneToOneRatio >= 0.70
    oneToOneCheck = oneToOneCheck and (mr.basesInOneToOne >= 5000)
    pvalCheck = mr.pval <= 0.01
    return (lengthCheck and missRateCheck and unAlignedCheck and oneToOneCheck and pvalCheck)

########################################################################
# parse matches for individual chromosomes
# each file in fileList represents the soma .all_match file for a chromosome
# Parse the matches, filter out the "bad matches", then from remaining matches
# take the best match for each contig for each chromosome
def parseChromosomeMatches(fileList):
    chromosomes = [getBaseName(f) for f in fileList]
    numChromosomes = len(chromosomes)
    filteredMatchLists = [parseMatchFileQuality(f) for f in fileList]
    return filteredMatchLists

########################################################################
def writeInfoFile(matchList, infoFileName):
    # For each quality match, write to infoFile:
    # contigName mapStart mapEnd forward NumSites misses missRate chiSquare medianSD lengthRatio
    infoFile = open(infoFileName, 'w')
    header = '%15s\t%15s\t%5s\t%5s\t%5s\t'%('contigName','chrom', 'start', 'end', 'Forwd')
    header += '%5s\t%10s\t%4s\t%4s\t%5s\t%8s\t%20s\t%10s\t%6s\t'%('sites','size','Cmiss','Omiss','missr','fragMiss','unalignedBaseRatio','score','pval')
    header += '%6s\t%6s\t%6s\t%6s\t%6s\t%10s\t%6s\n'%('meanSD','medSD','minSD','maxSD','lRat','1to2b', '1to1r')

    infoFile.write(header)
    for mr in matchList:
        line = '%15s\t%15s\t%5i\t%5i\t%5i\t'%(mr.contigId, mr.chromosome, mr.opStart, mr.opEnd, mr.forward)
        line += '%5i\t%10i\t%4i\t%4i\t%5.2f\t%8i\t%20.3f\t%10.2f\t%6.4f\t'%(mr.numContigSites, mr.size, mr.contigMisses, mr.opticalMisses, 100*mr.totalMissRate, mr.contigUnalignedFrags, mr.contigUnalignedBaseRatio, mr.score, mr.pval)
        line += '%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%10i\t%6.3f\n'%(mr.meanAlignedSD, mr.medianAlignedSD, mr.minAlignedSD, mr.maxAlignedSD, mr.lengthRatio, mr.basesInOneToOne, mr.basesInOneToOneRatio)
        infoFile.write(line)
    infoFile.close()

def writeInfoFile2(matchList, infoFileName):

    if not matchList:
        return

    fout = open(infoFileName, 'w')

    formatDict = {}
    formatDict[float] = lambda v: '%6.3f'%v

    fields = ['contigId', 
              'chromosome',
              'cStartIndex',
              'cEndIndex',
              'cStartBp',
              'cEndBp',
              'cAlignedBases',
              'opStartIndex',
              'opEndIndex',
              'opStartBp',
              'opEndBp',
              'opAlignedBases',
              'score']

    mr = matchList[0]
    fieldTypes = [type(mr.__getattribute__(f)) for f in fields]
    formatters = [formatDict.get(t, str) for t in fieldTypes]

    # Write Header
    fout.write('#' + ','.join(fields) + '\n')
    fd = formatDict
    for mr in matchList:
        fieldStrs = [formatter(mr.__getattribute__(f)) for f,formatter in zip(fields, formatters)]
        fout.write(','.join(fieldStrs) + '\n')
    fout.close()



########################################################################
# open a .all_match file, create a matchList of only quality matches,
# pickle the matchList to a file
def parseAndPickleQuality(matchFile, pickleFileName):
    matchList = parseMatchFileXMLQuality(matchFile)
    pickleFile = open(pickleFileName, 'w')
    cPickle.dump(matchList, pickleFile)
    pickleFile.close()

########################################################################
def parseAndPickle(matchFile, pickleFileName):
    matchList = parseMatchFileXML(matchFile)
    pickleFile = open(pickleFileName, 'w')
    cPickle.dump(matchList, pickleFile)
    pickleFile.close()

########################################################################
# Return -1 if mr1 < mr2 (meaning mr1 ranked higher)
# Return +1 if mr1 > mr2
# Return 0 if mr1 == mr2
# Compare basd on score, then...
# Compare based on number of misses,
# then median SD, then ChSquared, then length ratio
def cmpMatchResultsByQual(mr1, mr2):
    if mr1.score > mr2.score:
        return -1 # M1 ranked higher
    elif mr1.score < mr2.score:
        return 1 # M1 ranked lower
    elif mr1.contigUnalignedBaseRatio < mr2.contigUnalignedBaseRatio:
        return -1 #M1 ranked higher
    elif mr1.contigUnalignedBaseRatio > mr2.contigUnalignedBaseRatio:
        return 1 #M1 ranked lower
    elif mr1.totalMisses < mr2.totalMisses:
        return -1 #M1 ranked higher
    elif mr1.totalMisses > mr2.totalMisses:
        return 1
    elif mr1.medianAlignedSD < mr2.medianAlignedSD:
        return -1 #M1 ranked higher
    elif mr1.medianAlignedSD > mr2.medianAlignedSD:
        return 1
    elif mr1.lengthRatio > mr2.lengthRatio:
        return -1
    elif mr1.lengthRatio < mr2.lengthRatio:
        return 1
    else:
        return 0

########################################################################
# Return -1 if mr1 < mr2 (meaning mr1 is located on a lesser numbered chromosome or earlier in same chromosome)
# Return +1 if mr1 > mr2
# Return 0 if mr1 == mr2
# Compare based on number of misses,
# then median SD, then ChSquared, then length ratio
def cmpMatchResultsByChrom(mr1, mr2):
    p = re.compile(r'chr(\d+)')
    def getChromNum(chromString):
        mr = p.search(chromString)
        if not mr:
            return -1
        else:
            return int(mr.group(1))
    chrom1Num = getChromNum(mr1.chromosome)
    chrom2Num = getChromNum(mr2.chromosome)
    if chrom1Num < chrom2Num:
        return -1
    elif chrom1Num > chrom2Num:
        return 1
    elif mr1.opStart < mr2.opStart:
        return -1
    elif mr1.opStart > mr2.opStart:
        return 1
    elif mr1.opEnd < mr2.opEnd:
        return -1
    elif mr1.opEnd > mr2.opEnd:
        return 1
    else:
        return 0

########################################################################
def cmpMatchResultsBy1To1Ratio(mr1, mr2):
    if mr1.basesInOneToOneRatio > mr2.basesInOneToOneRatio:
        return -1
    elif mr1.basesInOneToOneRatio < mr2.basesInOneToOneRatio:
        return 1
    else:
        return 0

########################################################################
# Given a matchList:
# - Remove matches for contigs with less than 3 sites
# - Choose a set of non-overlapping alignments for each contig in the matchlist
# Return the refined matchList
# Rank by cmpMatchResultsByQual
# Added ranked contigs to the set in order, provided there is no overlap with a previously added contig
# NOTE: refining may be necessary becuase there are a lot of nearly identical quality matches to the same location, with
# a small difference in the number of misses or gaps, etc.
def refineMatchList(matchList):
    siteCutoff = 2
    # Create a dictionary from contigName to list of matches
    contigMatchDict = {}
    for mr in matchList:
        if mr.numContigSites < siteCutoff:
            continue
        if mr.contigId not in contigMatchDict: contigMatchDict[mr.contigId] = []
        contigMatchDict[mr.contigId].append(mr)
    refinedMatchList = []
    
    for contigId in sorted(contigMatchDict.keys()):
        contigMatches = contigMatchDict[contigId]
        if len(contigMatches) == 1:
            refinedMatchList.append(contigMatches[0])
            continue
        # This contig has multiple matches. Must resolve conflicts
        contigRefinedList = []
        sortedMatchList = sorted(contigMatches, cmp=cmpMatchResultsByQual)
        # Added match results to contigRefinedList sequentially,
        # provided that match result does not overlap with a result already
        # in the contigRefinedList
        for mr1 in sortedMatchList:
            hasOverlap = False
            for mr2 in contigRefinedList:
                # Check if mr1 and mr2 overlap
                if mr1.chromosome != mr2.chromosome:
                    continue
                if mr1.forward != mr2.forward:
                    continue
                # Being here means that mr1 and mr2 are on same chromosome
                # and in same orientation
                if not (mr1.opStart >= mr2.opEnd or mr1.opEnd <= mr2.opStart):
                    hasOverlap = True
                    break
            if not hasOverlap:
                contigRefinedList.append(mr1)
        refinedMatchList.extend(contigRefinedList) 

    return refinedMatchList

########################################################################
# Read XML file, generate list of MatchResult's, filter the list,
# and write results to .info and .pickle files.
def parse(xmlFileName):
    bn = getBaseName(xmlFileName)
    
    # parse the xmlFile for MatchResults
    # and pickle the list of matchResults
    pickleFileName = bn + '.pickle'
    parseAndPickle(xmlFileName, pickleFileName)

    # unpickle the matchList
    pickleFile = open(pickleFileName)
    ml = cPickle.load(pickleFile)
    pickleFile.close()

    # Write a human readable info files
    writeInfoFile(ml, bn+'.info')
    ml = sorted(ml, cmp=cmpMatchResultsByChrom)
    writeInfoFile(ml, bn+'.chromsort.info')

    pickleFileName = bn + '.filteredMatches.pickle'
    infoFileName = bn + '.filteredMatches.info'

    # Find high quality alignments from SOMA for each chromosome
    # Pickle the match list to a file
    parseAndPickleQuality(xmlFileName, pickleFileName)

    # Unpickle the quality match list
    pickleFile = open(pickleFileName)
    ml = cPickle.load(pickleFile)
    print 'Read %i matches from file %s'%(len(ml), pickleFileName)
    pickleFile.close()

    # Refine the matchLists, removing near duplicate matches for each contig
    rml = refineMatchList(ml)
    print 'Refined ml from %i matches to %i matches'%(len(ml),len(rml))

    # Organize output by contig
    contigToMatchList = {}
    finalMatchList = [] # Contigs with one or more quality alignment
    finalMatchListUnique = [] # Contigs that align to a unique location
    finalMatchListBest = [] # Best of the unique alignments
    for mr in rml:
        if mr.contigId not in contigToMatchList: contigToMatchList[mr.contigId] = []
        contigToMatchList[mr.contigId].append(mr)
    for contigId, contigMatchList in contigToMatchList.iteritems():
        finalMatchList.extend(contigMatchList)
        # If contig has only 1 alignment, append the match result to the unique match list
        if len(contigMatchList) == 1: finalMatchListUnique.append(contigMatchList[0])

    finalMatchListBest = [mr for mr in finalMatchListUnique if isQualityMatchA(mr)]

    print '%i contigs with significant alignments'%(len(contigToMatchList))
    print '%i contigs with only one significant alignment'%(len(finalMatchListUnique))
    print '%i contigs considered to be the best'%(len(finalMatchListBest))

    # Write the finalMatchListUnique to a pickle file
    pickleFile = open(bn + '.refinedUniqueMatchList.pickle', 'w')
    cPickle.dump(finalMatchListUnique, pickleFile)
    pickleFile.close()

    # Write the finalMatchList to a pickle file
    pickleFile = open(bn + '.refinedMatchList.pickle', 'w')
    cPickle.dump(finalMatchList, pickleFile)
    pickleFile.close()

    pickleFile = open(bn + '.refinedUniqueMatchListBest.pickle', 'w')
    cPickle.dump(finalMatchListBest, pickleFile)
    pickleFile.close()

    # Write the match lists to an info file
    finalMatchList = sorted(finalMatchList, cmp=cmpMatchResultsBy1To1Ratio)
    finalMatchListChromSort = sorted(finalMatchList, cmp=cmpMatchResultsByChrom)
    finalMatchListUniqueChromSort = sorted(finalMatchListUnique, cmp=cmpMatchResultsByChrom)
    finalMatchListBestChromSort = sorted(finalMatchListBest, cmp=cmpMatchResultsByChrom)
    writeInfoFile(finalMatchList, bn+'.refined.info')
    writeInfoFile(finalMatchListChromSort, bn+'.refined.chromsort.info')
    writeInfoFile(finalMatchListUnique, bn+'.refined.unique.info')
    writeInfoFile(finalMatchListUniqueChromSort, bn+'.refined.unique.chromsort.info')
    writeInfoFile(finalMatchListBestChromSort, bn+'.refined.unique.best.chromsort.info')

########################################################################
def plotFragLengthData():
    uniquePickleFile = open('refinedUniqueMatchList.pickle')
    ml = cPickle.load(uniquePickleFile)
    uniquePickleFile.close()
    # For the unique matches, plot the alignedOpticalFragLengths vs the alignedContigFragLengths (only the inner frags)
    alignedOpticalFragLengths = np.array([])
    alignedContigFragLengths = np.array([])
    alignedSD = np.array([])
    for mr in ml:
        alignedOpticalFragLengths = np.append(alignedOpticalFragLengths, mr.alignedOpticalLengths[1:-1])
        alignedContigFragLengths = np.append(alignedContigFragLengths, mr.alignedContigLengths[1:-1])
        alignedSD = np.append(alignedSD, np.sqrt(mr.alignedOpticalVariance[1:-1]))
    maxVal = np.max([np.max(alignedOpticalFragLengths), np.max(alignedContigFragLengths)])
    plt.errorbar(alignedContigFragLengths, alignedOpticalFragLengths,yerr=alignedSD,linestyle='None',marker='.')
    plt.plot([0, maxVal], [0, maxVal],'g')
    plt.axis('equal')
    plt.show()

#############################################
# Print alignments in the match list out to file (or stdout)
def printAlignments(ml, fout = sys.stdout):
    for mr in ml:
        fout.write('*'*50+'\n')
        mr.printAlignment(fout)

########################################################################
if __name__ == '__main__':
    xmlFileName = sys.argv[1]
    parse(xmlFileName)
