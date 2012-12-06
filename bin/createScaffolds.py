# Create scaffolds from unique soma alignments
import importSomaData
import glob
import numpy as np


# Return a tiling of the matchList to
# the chromosomes. Sort the matchList
# by size. Place a contig to a chromosome if there
# is not overlap with a previously place contig
# The tiling is simply represented as a dictionary:
# keys: chromosome, value: list of sorted MatchResults for chromosome
def createTiling(matchList, allowOverlaps=False):
    # Sort the matchList by size
    matchList = sorted(matchList, key=lambda mr: mr.size)[::-1]
    chromToML = {}
    noOverlap = lambda mr1, mr2: (mr1.opStart >= mr2.opEnd) or (mr1.opEnd <= mr2.opStart)
    for mr in matchList:
        if mr.chromosomeName not in chromToML:
            chromToML[mr.chromosomeName] = [mr]
        else:
            numOverlaps = sum(int(not noOverlap(mr, mr2)) for mr2 in chromToML[mr.chromosomeName])
            if (numOverlaps==0 or allowOverlaps):
                chromToML[mr.chromosomeName].append(mr)
    for chrom in chromToML.keys():
        chromToML[chrom] = sorted(chromToML[chrom], key = lambda mr: mr.opStart)
    return chromToML


def writeScaffoldsFile(chromToML, chromStats, overallStats, scaffFileName):
    fout = open(scaffFileName, 'w')

    # Write overall scaffold stats
    fout.write('Total Optical Map Size: %i\n'%overallStats['opMapSize'])
    fout.write('Contigs Placed in Scaffold: %i\n'%overallStats['contigsInScaff'])
    fout.write('Covered bases: %i\n'%overallStats['coveredSize'])
    fout.write('Covered fraction: %8.5f\n'%overallStats['coveredRatio'])
    fout.write('*'*50 + '\n\n')

    chromList = sorted(chromToML.keys())
    for chrom in chromList:
        stats = chromStats[chrom]
        header = '*'*50 + '\n' + chrom + '\n' + '*'*50 + '\n'
        header += 'chromSize: %i coveredSize: %i coveredRatio: %f\n'%(stats['size'], stats['coveredSize'], stats['coveredRatio'])
        header += '%20s\t%10s\t%10s\t%10s\t%10s'%('contigId','length','start','end', 'gap after') + '\n'
        fout.write(header)
        ml = chromToML[chrom]
        gaps = [ml[i].opStartBp - ml[i-1].opEndBp for i in range(1,len(ml))]
        gaps.append(0)
        fout.write('\n'.join(['%20s\t%10i\t%10i\t%10i\t%10i'%(ml[i].contigId,ml[i].size,ml[i].opStartBp, ml[i].opEndBp, gaps[i]) for i in range(len(ml))]))
        fout.write('\n')
    fout.close()

def createScaffolds(matchListPickle, opticalMapFiles, outFile, allowOverlaps=False):

    # Read matches
    contigToMatchList = importSomaData.getContigToMatchList(matchListPickle)
    ml = []
    for contigId, contigMatches in contigToMatchList.iteritems():
        ml.extend(contigMatches)
    print 'Read %i matches from pickle files %s'%(len(ml), matchListPickle)

    # Read optical maps
    opticalMapDict = {}
    for oFile in opticalMapFiles:
        opticalMapDict[oFile] = importSomaData.readOpticalMap(oFile)

    # Create a tiling of the Match Results
    chromToML = createTiling(ml, allowOverlaps=allowOverlaps)

    overallStats = {}

    # For each tiling, computing the starting and ending location (in bp) of each matchResult with 
    # respect to the optical map.
    chromStats = {}
    opMapSize = np.sum([np.sum(om[:,0]) for om in opticalMapDict.values()]) 
    totalCoveredSize = 0
    totalContigsInScaff = 0
    for oMap in chromToML.keys():
        opticalFragLengths = opticalMapDict[oMap][:,0]
        cumulativeSum = np.cumsum(opticalFragLengths)
        fragStart = np.concatenate(([0], cumulativeSum[0:-1]))
        ml = chromToML[oMap]
        for mr in ml:
            alignedContigFragLengths = mr.alignedContigLengths
            mr.opStartBp = fragStart[mr.opStart+1]-alignedContigFragLengths[0] 
            mr.opEndBp = fragStart[mr.opEnd]+alignedContigFragLengths[-1]
        totalContigsInScaff += len(ml)
        chromSize = np.sum(opticalFragLengths)
        coveredSize = np.sum([mr.size for mr in ml])
        totalCoveredSize += coveredSize
        coveredRatio = coveredSize/float(chromSize)
        # Note: The coveredSize and coveredRatio stats will be skewed if overlaps are allowed
        chromStats[oMap] = dict([('size',chromSize), ('coveredSize', coveredSize), ('coveredRatio',coveredRatio)])
    overallStats['opMapSize'] = opMapSize
    overallStats['contigsInScaff'] = totalContigsInScaff
    overallStats['coveredSize'] = totalCoveredSize
    overallStats['coveredRatio'] =  float(totalCoveredSize)/opMapSize
    writeScaffoldsFile(chromToML, chromStats, overallStats, outFile)

def main():
    pass

if __name__=='__main__':
    main()
