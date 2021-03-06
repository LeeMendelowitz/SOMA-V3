# Summarize the alignment status for each contig
# Want to quickly display information such as:
# contig, contigSize, number restriction sites,
# aligned uniquely, numberOfQualityalignments, number of chromsomes aligned to
import importSomaData
import numpy as np
import SOMAMap
import sys
import os

def getBaseName(inFile):
    (dirName, fileName) = os.path.split(inFile)
    (fileBaseName, fileExtension)=os.path.splitext(fileName)
    return fileBaseName

# Get contigs that have a significantly large fragment,
# making these contigs particularly unique and (hopefully) alignable
'''
def getSpecialContigs():
    silicoFile = 'contigs.silico'
    contigDict = SOMAMap.readMaps(silicoFile)
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
'''

# Summarize alignment status for the contigs in contigList
# If no contigs are provided, summarize results for all contigs 
def summarizeContigStatus(outputPfx, contigMatchDict, contigMapDict, contigList = None):

    outputFile = '%s.contigAlignmentSummary.csv' %outputPfx
    uniqueMatches = set([c for c,matchlist in contigMatchDict.iteritems() if len(matchlist)==1])

    # write to a summary file
    fout = open(outputFile, 'w')
    header = 'contig,size,numFrags,numAlignments,alignedUniquely,numChroms,chroms'
    fout.write(header + '\n')
    if not contigList:
        contigList = contigMapDict.keys()
    # sort contigs by size descending
    contigList = sorted(contigList, key = lambda s: contigMapDict[s].length)[::-1]
    for c in contigList:
        contigMap = contigMapDict[c]
        numAlignments = len(contigMatchDict[c]) if c in contigMatchDict else 0
        alignedUniquely = c in uniqueMatches
        numChroms = 0
        chroms = ''

        if numAlignments:
            chromSet = set([mr.chromosome for mr in contigMatchDict[c]])
            numChroms = len(chromSet)
            chroms = ';'.join(chromSet)
        fields   = [c, 
                    contigMap.length,
                    contigMap.numFrags,
                    numAlignments,
                    alignedUniquely,
                    numChroms,
                    chroms
                    ]
        lineOut = ','.join(str(f) for f in fields)
        fout.write(lineOut + '\n')
    fout.close()

def readContigSummary(summaryFile):
    raise Exception('TO BE IMPLEMENTED!')
#    fin = open(summaryFile)
#    fieldNames = ['contig','size','numSites','numAlignments','alignedUniquely','alignedBest','numChroms','chroms']
#    def readChroms(s):
#        if s:
#            return s.split(';')
#        return []
#    str2bool = lambda s: s=='True'
#    types = [str, int, int, int,str2bool,str2bool,int,readChroms]
#    recDict = {}
#    fin.next()
#    fieldGen = (l.strip().split(',') for l in fin)
#    for fields in fieldGen:
#        rec = dict((fn, t(f)) for fn,t,f in zip(fieldNames,types,fields))
#        recDict[rec['contig']] = rec
#
#    return recDict

def main():
    bn = getBaseName(sys.argv[1])
    silicoFile = sys.argv[2]
    summarizeContigStatus(bn, silicoFile)


if __name__ == '__main__':
    main()
