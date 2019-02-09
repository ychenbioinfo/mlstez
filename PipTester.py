#!/usr/bin/env python
import sys
import Parameters
import ProjectEnviroment
import SeqAlignParallel
from time import time

def FileModule():
    parameter = Parameters.Parameters()
    parameter.openfile("config.ini")
    parameter.update()
    projEnv = ProjectEnviroment.ProjectEnviroment(parameter)
    projEnv.loadFiles()
    barcodes = projEnv.Barcodes
    primers = projEnv.Primers
    seqs = projEnv.Seqs
    t0 = time()
    projEnv.alignSeqs()
    alignseqs = projEnv.AlignedSeqs
    print ("Sort sequencing...")
    isokay, error = projEnv.sortAlignedSeqs()
    if(not isokay):
        print(("Error: %s" %error))
    sortedseqs = projEnv.SortedSeqs
    print ("Calculate length...")
    isokay, error = projEnv.locusLengths()
    if(not isokay):
        print(("Error: %s" %error))
    locusRange = projEnv.locusLengthRange
    print ("Generate Consensus Sequences...")
    isokay, error = projEnv.GenerateConsensus()
    if(not isokay):
        print(("Error: %s" %error))
    consSeqs = projEnv.consSeqs
    projEnv.DumpUnmappedReads()
    t1 = time()
    print(("Time spends %.2fs." %(t1-t0)))


if(__name__ == '__main__'):
    FileModule()
    