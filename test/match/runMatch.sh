#!/bin/bash -x
BIN_DIR=../../bin
MATCH_BIN_DEBUG=$BIN_DIR/debug/match
MATCH_BIN_RELEASE=$BIN_DIR/release/match

# TEST RELEASE
$MATCH_BIN_RELEASE --output test100WithError --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 1000 --numThreads 4 contigs.silico map100.opt
#$MATCH_BIN_RELEASE --output test100_maxMatch1 --maxMatchesPerContig=1 --sdMin=4 --sdMax=4 --numPermutationTrials 0 --numThreads 4 contigs.silico map100.opt
$MATCH_BIN_RELEASE --output test100_maxMatch1_noBoundary --maxMatchesPerContig=1  --noBoundaries --sdMin=4 --sdMax=4 --numPermutationTrials 0 --numThreads 4 contigs.silico map100.opt

# TEST DEBUG
#$MATCH_BIN_DEBUG --output test100WithError --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 0 --numThreads 1 contigs.silico map100.opt
#$MATCH_BIN_DEBUG --output test100WithError --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 0 --numThreads 1 contig2.silico map100.opt
#$MATCH_BIN_DEBUG --output contig14_noBoundary --maxMatchesPerContig=1  --noBoundaries --sdMin=4 --sdMax=4 --numPermutationTrials 0 --numThreads 1 contig14.silico map100.opt
