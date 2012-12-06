#!/usr/bin/env python
import sys
import subprocess
import re
import os
from Bio import SeqIO
import argparse

# Usage: make_silico.py fasta enzyme outfile
class Site:
    def __init__(self,loc):
       self.loc = loc

    # Return string representation
    def str(self):
        return str(self.loc)

# find restriction sites 
# dna is a contig
def findSites(dna, recSeq):
    sites = []
    lr = len(recSeq)
    for loc in range(0, len(dna)):
        if(dna[loc:loc+lr] == recSeq):
            sites.append(Site(loc+1))
    return sites
   
def createInsilico(fout, dna, contigId, recSeq):
    ldna = len(dna)
    sites = findSites(dna, recSeq)
    numSites = len(sites)
    if (ldna > 0 and numSites > 0):
        outString = ' '.join([contigId,str(ldna), str(numSites)])+'\n'
        outString += ';'.join([site.str() for site in sites])
        fout.write(outString + '\n')

# Return arguments:
#   fasta, enzyme, outfile
def parseArgs():
    parser = argparse.ArgumentParser(description="Create in-silico optical map from FASTA")
    parser.add_argument('fasta')
    parser.add_argument('enzyme')
    parser.add_argument('outfile')
    args = parser.parse_args()
    return args

def getRecSeq(enzyme):
    repbase = 'repbase.staden'
    if not os.path.exists(repbase):
        raise RuntimeError('make_silico ERROR: Cannot find %s'%repbase)

    # Get the enzyme recognition pattern
    mLines = [l.strip() for l in open(repbase) if enzyme in l]
    if not mLines:
        raise RuntimeError('No Match for enzyme %s in %s\n'%(enzyme, repbase))
    if len(mLines)>1:
        raise RuntimeError('Found %i matches for enzyme %s in %s\n'%(len(mLines), enzyme, repbase))
    grepresult = mLines[0]
    name, recSeq = grepresult.split('/')[0:2]
    recSeq = recSeq.replace('\'','')

    # Check that the resctriction_seq is not degenerate
    if(re.search('[^acgtACGT]', recSeq)):
        raise RuntimeError('ERROR: %s has a degenerate recognition sequence %s'%(enzyme, recSeq))

    # Convert recognition sequence to upper case
    recSeq = recSeq.upper()
    return recSeq
   
def makeInsilicoMain(fasta, enzyme, outFile):
    recSeq = getRecSeq(enzyme)
    fin = open(fasta)
    fout = open(outFile, 'w')
    for rec in SeqIO.parse(fin, 'fasta'):
        dna = str(rec.seq).upper()
        contigId = rec.id
        createInsilico(fout, dna, contigId, recSeq)
    fin.close()
    fout.close()

def main():
    args = parseArgs()
    fasta = args.fasta
    enzyme = args.enzyme
    outFile = args.outfile
    makeInsilicoMain(fasta, enzyme, outFile)


if __name__ == '__main__':
    main()
