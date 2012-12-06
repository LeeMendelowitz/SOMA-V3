# Summarize the alignment status for each contig
# Want to quickly display information such as:
# contig, contigSize, number restriction sites,
# aligned uniquely, aligned best, numberOfQualityalignments, number of chromsomes aligned to
#import importSomaData
import numpy as np

import sys
import os

def getBaseName(inFile):
    (dirName, fileName) = os.path.split(inFile)
    (fileBaseName, fileExtension)=os.path.splitext(fileName)
    return fileBaseName

# Get contigs that have a significantly large fragment,
# making these contigs particularly unique and (hopefully) alignable
def getSpecialContigs():
    import importSomaData
    silicoFile = 'contigs.silico'
    contigDict = importSomaData.readSilicoFile(silicoFile)
    allFragLengths = np.concatenate([c.frag for c in contigDict.itervalues()])/1000.0
    allFragLengths = np.sort(allFragLengths)
    numFrags = allFragLengths.shape[0]
    innerFragLengths = np.concatenate([c.frag[1:-1] for c in contigDict.itervalues()])/1000.0
    innerFragLengths = np.sort(innerFragLengths)
    fracBelow = np.arange(1.0,numFrags+1)/numFrags
    findPercentile = lambda c: allFragLengths[fracBelow >= c][0]
    # Find and return contigs that have a fragment that is above the 0.995 percentile
    # in length
    cutoff = findPercentile(0.995)
    cList = []
    for c, contigData in contigDict.iteritems():
        if np.max(contigData.frag/1000.0) > cutoff:
            cList.append(c)
    return cList

# Summarize alignment status for the contigs in contigList
# If no contigs are provided, summarize results for all contigs 
def summarizeContigStatus(bn, silicoFile, contigList = None):
    import importSomaData

    outputFile = bn + '.contigAlignmentSummary.csv' 
    # Read contig id's, size, and number of sites from the .silico file
    contigDataDict = importSomaData.readSilicoFile(silicoFile)

    # Read soma alignments
    matchFile = bn + '.refinedMatchList.pickle'
    uniqueMatchFile = bn + '.refinedUniqueMatchList.pickle'
    bestMatchfile = bn + '.refinedUniqueMatchListBest.pickle'
    
    print 'Reading contigMatchDict....'
    contigMatchDict = importSomaData.getContigToMatchList(matchFile)
    print 'DONE.'
    uniqueMatches = set([c for c,matchlist in contigMatchDict.iteritems() if len(matchlist)==1])
    bestMatches = importSomaData.getContigToMatchList(bestMatchfile)

    # write to a summary file
    fout = open(outputFile, 'w')
    header = 'contig,size,numSites,numAlignments,alignedUniquely,alignedBest,numChroms,chroms'
    fout.write(header + '\n')
    if not contigList:
        contigList = contigDataDict.keys()
    # sort contigs by size descending
    contigList = sorted(contigList, key = lambda s: contigDataDict[s].length)[::-1]
    for c in contigList:
        contigData = contigDataDict[c]
        numAlignments = len(contigMatchDict[c]) if c in contigMatchDict else 0
        alignedUniquely = c in uniqueMatches
        alignedBest = c in bestMatches
        numChroms = 0
        chroms = ''

        if numAlignments:
            chromSet = set([mr.chromosomeName for mr in contigMatchDict[c]])
            numChroms = len(chromSet)
            chroms = ';'.join(chromSet)
        lineOut = ','.join(map(str,[c, contigData.length, contigData.numSites,
                  numAlignments, alignedUniquely, alignedBest, numChroms, chroms]))
        fout.write(lineOut + '\n')
    fout.close()

def readContigSummary(summaryFile):
    fin = open(summaryFile)
    fieldNames = ['contig','size','numSites','numAlignments','alignedUniquely','alignedBest','numChroms','chroms']
    def readChroms(s):
        if s:
            return s.split(';')
        return []
    str2bool = lambda s: s=='True'
    types = [str, int, int, int,str2bool,str2bool,int,readChroms]
    recDict = {}
    fin.next()
    fieldGen = (l.strip().split(',') for l in fin)
    for fields in fieldGen:
        rec = dict((fn, t(f)) for fn,t,f in zip(fieldNames,types,fields))
        recDict[rec['contig']] = rec

    return recDict



def main():
    bn = getBaseName(sys.argv[1])
    silicoFile = sys.argv[2]
    summarizeContigStatus(bn, silicoFile)



if __name__ == '__main__':
    main()
