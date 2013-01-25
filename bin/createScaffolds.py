# Create scaffolds from unique soma alignments
import numpy as np
import itertools


# Return a tiling of the matchList to
# the chromosomes. Sort the matchList
# by size. Place a contig to a chromosome if there
# is not overlap with a previously place contig
# The tiling is simply represented as a dictionary:
# keys: chromosome, value: list of sorted MatchResults for chromosome
def createTiling(matchList, allowOverlaps=False):
    # Sort the matchList by size
    matchList = sorted(matchList, key=lambda mr: mr.cAlignedBases)[::-1]
    chromToML = {}
    noOverlap = lambda mr1, mr2: (mr1.opStartBp <= mr2.opEndBp) or (mr1.opEndBp <= mr2.opStartBp)
    for mr in matchList:
        if mr.chromosome not in chromToML:
            chromToML[mr.chromosome] = [mr]
        else:
            if allowOverlaps:
                chromToML[mr.chromosome].append(mr)
                continue
            numOverlaps = sum(int(not noOverlap(mr, mr2)) for mr2 in chromToML[mr.chromosome])
            if (numOverlaps==0):
                chromToML[mr.chromosome].append(mr)
    for chrom in chromToML.keys():
        chromToML[chrom] = sorted(chromToML[chrom], key = lambda mr: mr.opStartBp)
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

    chromHeaderTemplate =  '*'*50 + '\n{chromId}\n' + '*'*50 + '\n'
    chromHeaderTemplate += 'chromSize: {stats[length]} coveredSize: {stats[coveredSize]} coveredRatio: {stats[coveredRatio]:6.3f}\n'

    # (name, header spec, record spec)
    fields = [ ('contigId', '{0:20s}', '{match.contigId:20s}' ),
               ('length', '{0:10s}', '{match.cAlignedBases:10d}'),
               ('start', '{0:10s}', '{match.opStartBp:10d}'),
               ('end', '{0:10s}', '{match.opEndBp:10d}'),
               ('gap after', '{0:10s}', '{gap:10d}')
             ]

    chromHeaderTemplate += '\t'.join(headerSpec.format(name) for name,headerSpec,recordSpec in fields) + '\n'
    recordTemplate = '\t'.join(recordSpec for name,headerSpec,recordSpec in fields) + '\n'

    for chrom in chromList:
        stats = chromStats[chrom]
        header = chromHeaderTemplate.format(chromId = chrom, stats=stats)
        #header = '*'*50 + '\n' + chrom + '\n' + '*'*50 + '\n'
        #header += 'chromSize: %i coveredSize: %i coveredRatio: %f\n'%(stats['size'], stats['coveredSize'], stats['coveredRatio'])
        #header += '%20s\t%10s\t%10s\t%10s\t%10s'%('contigId','length','start','end', 'gap after') + '\n'
        fout.write(header)
        ml = chromToML[chrom]
        gaps = [ml[i].opStartBp - ml[i-1].opEndBp for i in range(1,len(ml))] + [0]
        assert(len(ml) == len(gaps))
        for mr, gapAfter in itertools.izip(ml, gaps):
            record = recordTemplate.format(match = mr, gap=gapAfter)
            fout.write(record)
    fout.close()

#
# contigMatchDict: contigId to list of matches for that contig
def createScaffolds(contigMatchDict, opticalMapDict, outFile, allowOverlaps=False):

    ml = [mr for contigId, contigMatches in contigMatchDict.iteritems() for mr in contigMatches]
    # Create a tiling of the Match Results
    chromToML = createTiling(ml, allowOverlaps=allowOverlaps)

    overallStats = {}
    chromStats = {}

    opMapSize = sum(omap.length for omap in opticalMapDict.itervalues())
    totalCoveredSize = 0
    totalContigsInScaff = 0
    for omapId, chromML in chromToML.iteritems():
        omap = opticalMapDict[omapId]
        totalContigsInScaff += len(chromML)
        coveredSize = np.sum([mr.opEndBp - mr.opStartBp for mr in chromML])
        totalCoveredSize += coveredSize
        coveredRatio = coveredSize/float(omap.length)
        # Note: The coveredSize and coveredRatio stats will be skewed if overlaps are allowed
        chromStats[omapId] = dict([('length', omap.length), ('coveredSize', coveredSize), ('coveredRatio',coveredRatio)])
    overallStats['opMapSize'] = opMapSize
    overallStats['contigsInScaff'] = totalContigsInScaff
    overallStats['coveredSize'] = totalCoveredSize
    overallStats['coveredRatio'] =  float(totalCoveredSize)/opMapSize
    writeScaffoldsFile(chromToML, chromStats, overallStats, outFile)

def main():
    pass

if __name__=='__main__':
    main()
