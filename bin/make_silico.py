#!/usr/bin/env python
import sys
import subprocess
import re
import os
from Bio import SeqIO
import argparse
from rebase import rebase

# Usage: make_silico.py fasta enzyme outfile

# find restriction sites 
# dna is a contig
def findSites(dna, recSeq):
    sites = []
    lr = len(recSeq)
    sites = [loc for loc in xrange(len(dna)) if dna[loc:loc+lr]==recSeq]
    return sites
  
##############################################
# Write an insilico restriction map for the given contig
def writeInsilico(fout, dna, contigId, recSeq):
    ldna = len(dna)
    sites = findSites(dna, recSeq)
    numSites = len(sites)
    frags = []

    # compute restriction fragments from sites
    if sites:
        frags = [sites[0]]
        frags.extend(sites[i+1] - sites[i] for i in xrange(0, numSites-1))
        frags.append(ldna - sites[-1])

    outS = [contigId, str(ldna), str(len(frags))]
    outS.extend(str(frag) for frag in frags)
    outS = '\t'.join(outS)
    fout.write(outS + '\n')

# Return arguments:
#   fasta, enzyme, outfile
def parseArgs():
    parser = argparse.ArgumentParser(description="Create in-silico optical map from FASTA")
    parser.add_argument('fasta', help = 'Path to input FASTA file')
    parser.add_argument('enzyme', help = 'Enzyme name')
    parser.add_argument('outfile', help = 'Output file.')
    args = parser.parse_args()
    return args

def getRecSeq(enzyme):
   
    enzyme = enzyme.upper()
    # Get the enzyme recognition pattern
    if enzyme not in rebase:
        raise RuntimeError('No Match for enzyme %s in rebase'%enzyme)

    recSeqs = rebase[enzyme]

    if len(recSeqs)>1:
        raise RuntimeError('Found %i matches for enzyme %s\n'%(len(mLines), enzyme))

    recSeq = recSeqs[0]

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
        writeInsilico(fout, dna, contigId, recSeq)
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
