from MUSCython import MultiStringBWTCython as msbwt, LCPGen
from time import clock
import logging
from multiprocessing import Pool
from itertools import izip, repeat, izip_longest
from random import choice
from os import remove
from string import maketrans
import cython
import numpy as np
cimport numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plot

cdef int K = 25
cdef list BASES = ['A', 'C', 'G', 'T']
cdef int READ_LEN = 100
cdef int THRESH = 1
cdef int HI_THRESH = 5
cdef int NUM_PROCS = 4
# PATH = '/playpen/sgreens/ecoli/'
# BWT_PATH = '/playpen/sgreens/ecoli/bwt/'
PATH = '/playpen/sgreens/ecoli/msbwt20/'
BWT_PATH = '/playpen/sgreens/ecoli/msbwt20/rle_bwt2/'
# PATH = '/playpen/sgreens/ecoli/msbwtpaired20/'
# BWT_PATH = '/playpen/sgreens/ecoli/msbwtpaired20/bwt/'
# PATH = '/playpen/sgreens/ecoli/msbwtNoN20/'
# BWT_PATH = '/playpen/sgreens/ecoli/msbwtNoN20/bwt/'
# PATH = '/playpen/sgreens/fake/'
# BWT_PATH = '/playpen/sgreens/fake/bwt/'

BASE_NUMS = {'A': 1, 'C': 2, 'G': 3, 'T': 5}
NUM_TO_BASE = ['A', 'C', 'G', 'N', 'T', '$']

TRAN_TAB = maketrans('ACGT', 'TGCA')

def reverseComplement(seq):
    return seq[::-1].translate(TRAN_TAB)


def fastaParser(fasta, int start, int end):
    cdef long _
    with open(fasta, 'r') as fp:
        for _ in xrange(start*4):
                fp.next()
        for _ in xrange(start, end):
            yield fp.next(), fp.next(), fp.next(), fp.next()


def superCorrect(wn_nw):
    begin = clock()
    workerNum = wn_nw[0]
    numWorkers = wn_nw[1]
    bwt = msbwt.loadBWT(BWT_PATH, False)
    cdef np.ndarray[np.uint8_t] trusted = np.empty(READ_LEN-K+1, dtype=np.uint8)
    cdef np.ndarray[np.uint8_t] corrected = np.empty(READ_LEN, dtype=np.uint8)
    cdef np.ndarray[np.uint64_t] counts = np.empty(READ_LEN-K+1, dtype=np.uint64)
    cdef np.ndarray[np.uint64_t] revCounts = np.empty(READ_LEN-K+1, dtype=np.uint64)
    cdef np.ndarray[np.uint64_t] newCounts = np.empty(READ_LEN-K+1, dtype=np.uint64)
    cdef np.ndarray[np.int64_t] _ = np.empty(READ_LEN-K+1, dtype=np.int64)
    cdef int readsPerWorker = (bwt.getSymbolCount(0) / numWorkers) + 1  # reads per worker per file
    cdef long lo, hi
    cdef bytes readName, origRead, plus, qual, base, bestBase
    cdef bint changeMade
    print 'Worker' , workerNum, '/', numWorkers, 'doing', readsPerWorker*workerNum, 'through', min(readsPerWorker * (workerNum+1), bwt.getSymbolCount(0))
    fileSuffix = ''
    cdef int readsDone = 0
    with open(PATH + 'corrected' + fileSuffix + '_' + str(workerNum) + '.fastq', 'w') as fp:
        for readName, origRead, plus, qual in \
                fastaParser('/playpen/sgreens/ecoli/EAS20_8/cov20' + fileSuffix + '.txt',
                            readsPerWorker*workerNum,
                            min(readsPerWorker * (workerNum+1), bwt.getSymbolCount(0))):
            readsDone += 1
            # if readName == '@EAS20_8_6_1_51_932/1\n' or readName == '@EAS20_8_6_1_51_932/2\n':
            #     print origRead
            #     return
            if not readsDone % 1000:
                print 'Worker', workerNum, 'working on read', readsPerWorker*workerNum + readsDone, '/', bwt.getSymbolCount(0)
                print origRead
            read = list(origRead[:-1])
            revCounts, _ = bwt.countStrandedSeqMatches(reverseComplement(origRead[:-1]), K)
            revCounts = np.flipud(revCounts)
            trusted = (revCounts > 0).astype(np.uint8)
            corrected.fill(0)
            if False in trusted:
                counts, _ = bwt.countStrandedSeqMatches(origRead[:-1], K)
                trusted |= counts > HI_THRESH
                changeMade = True
                while changeMade:
                    changeMade = False
                    if False in trusted:
                        for i in xrange(len(trusted)-1):
                            if trusted[i] or read[i+K] == 'N':
                                if not trusted[i+1] or read[i+K] == 'N':  # err at read[i+K]
                                    bestSupport = 0
                                    for base in BASES:
                                        lo, hi = bwt.findIndicesOfStr(''.join(read[i+1:i+K]) + base)
                                        rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[i+1:i+K]) + base))
                                        if hi - lo + rcHi - rcLo > bestSupport:
                                            bestSupport = hi - lo + rcHi - rcLo
                                            bestBase = base
                                    if bestBase != read[i+K]:
                                        changeMade = True
                                        corrected[i+K] = True
                                        read[i+K] = bestBase
                                    if corrected[i+K]:
                                        kmersEnd = min(len(read), i+2*K)
                                        newCounts, _ = bwt.countStrandedSeqMatches(''.join(read[i+1:kmersEnd]), K)
                                        trusted[i+1:kmersEnd-K+1] = newCounts > HI_THRESH
                                        counts[i+1:kmersEnd-K+1] = newCounts
                                        if False in trusted[i+1:kmersEnd-K+1]:
                                            newCounts, _ = bwt.countStrandedSeqMatches(
                                                reverseComplement(''.join(read[i+1:kmersEnd])), K)
                                            newCounts = np.flipud(newCounts)
                                            trusted[i+1:kmersEnd-K+1] |= newCounts > 0
                            elif trusted[i+1] and not corrected[i]:  # err at read[i]
                                bestSupport = 0
                                for base in BASES:
                                    lo, hi = bwt.findIndicesOfStr(base + ''.join(read[i+1:i+K]))
                                    rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(base + ''.join(read[i+1:i+K])))
                                    if hi - lo + rcHi - rcLo > bestSupport:
                                        bestSupport = hi - lo + rcHi - rcLo
                                        bestBase = base
                                if bestBase != read[i]:
                                    changeMade = True
                                    corrected[i] = True
                                    read[i] = bestBase
                                if corrected[i]:
                                    kmersStart = max(0, i-K)
                                    newCounts, _ = bwt.countStrandedSeqMatches(''.join(read[kmersStart:i+K]), K)
                                    trusted[kmersStart:i+1] = newCounts > HI_THRESH
                                    counts[kmersStart:i+1] = newCounts
                                    if False in trusted[kmersStart:i+1]:
                                        newCounts, _ = bwt.countStrandedSeqMatches(
                                            reverseComplement(''.join(read[kmersStart:i+K])), K)
                                        newCounts = np.flipud(newCounts)
                                        trusted[kmersStart:i+1] |= newCounts > 0
            fp.write(readName)
            fp.write(''.join(read) + '\n')
            fp.write(plus)
            fp.write(qual)
    print 'Worker', workerNum, 'finished in', clock() - begin, 'seconds'


def main():
    NUM_WORKERS = 12
    pool = Pool(NUM_PROCS)
    pool.map(superCorrect, izip(range(NUM_WORKERS), repeat(NUM_WORKERS)))
    print 'All workers finished'
    fileSuffix = ''
    with open(PATH + 'corrected' + fileSuffix + '.fastq', 'w') as outFile:
        for i in xrange(NUM_WORKERS):
            fileName = PATH + 'corrected' + fileSuffix + '_' + str(i) + '.fastq'
            with open(fileName, 'r') as inFile:
                for line in inFile:
                    outFile.write(line)
            remove(fileName)


if __name__ == '__main__':
    main()