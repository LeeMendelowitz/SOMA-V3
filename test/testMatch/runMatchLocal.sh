#!/bin/bash
#/cbcb/project-scratch/lmendelo/soma/test/testMatch
BIN_DIR=../../bin
MATCH_BIN=$BIN_DIR/match_d
PTRIALS=0
THREADS=1
$MATCH_BIN --output test100WithError --local --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials $PTRIALS --numThreads $THREADS 100_WithError.silico 100.opt
$MATCH_BIN --output test100NoError --local --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials $PTRIALS --numThreads $THREADS 100_NoError.silico 100.opt
