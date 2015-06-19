from MUSCython import MultiStringBWTCython as msbwt, LCPGen
from time import clock
from tempfile import NamedTemporaryFile
import logging
from multiprocessing import Pool
from os import remove
from string import maketrans
from functools import partial
import cython
import numpy as np
cimport numpy as np

cdef list BASES = ['A', 'C', 'G', 'T']

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


def correct(inFastqFile, bwtDir, k, hiThresh, numProcesses, processNum):
    begin = clock()
    bwt = msbwt.loadBWT(bwtDir, False)
    cdef np.ndarray[np.uint8_t] corrected
    cdef np.ndarray[np.uint8_t] trusted
    cdef np.ndarray[np.uint64_t] counts
    cdef np.ndarray[np.uint64_t] revCounts
    cdef np.ndarray[np.uint64_t] newCounts
    cdef np.ndarray[np.int64_t] _
    cdef int readsPerProcess = (bwt.getSymbolCount(0) / numProcesses) + 1
    cdef int kmersEnd, kmersStart
    cdef long lo, hi, rcLo, rcHi, bestSupport
    cdef bytes readName, origRead, plus, qual, base, bestBase
    cdef bint changeMade
    tmpFile = NamedTemporaryFile(suffix='.fastq', delete=False)
    logging.info('Process %d / %d: correcting reads %d through %d', processNum, numProcesses,
                 readsPerProcess*processNum, min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))-1)
    cdef long readsDone = 0
    with tmpFile:
        for readName, origRead, plus, qual in \
                fastaParser(inFastqFile, readsPerProcess*processNum,
                            min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))):
            read = list(origRead[:-1])
            revCounts, _ = bwt.countStrandedSeqMatches(reverseComplement(origRead[:-1]), k)
            revCounts = np.flipud(revCounts)
            trusted = (revCounts > 0).astype(np.uint8)
            corrected = np.zeros(len(read), dtype=np.uint8)
            readsDone += 1
            if not readsDone & 1023:
                logging.debug('Finished %d reads', readsDone)
            if False in trusted:
                counts, _ = bwt.countStrandedSeqMatches(origRead[:-1], k)
                trusted |= counts > hiThresh
                changeMade = True
                while changeMade:
                    changeMade = False
                    if False in trusted:
                        for i in xrange(len(trusted)-1):
                            if trusted[i] or read[i+k] == 'N':
                                if not trusted[i+1] or read[i+k] == 'N':  # err at read[i+k]
                                    bestSupport = 0
                                    newLo, newHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[i+1:i+k])))
                                    for base in BASES:
                                        rcLo = bwt.getOccurrenceOfCharAtIndex(BASE_NUMS[base.translate(TRAN_TAB)], newLo)
                                        rcHi = bwt.getOccurrenceOfCharAtIndex(BASE_NUMS[base.translate(TRAN_TAB)], newHi)
                                        # rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[i+1:i+k]) + base))
                                        # if hi - lo + rcHi - rcLo > bestSupport:
                                        if rcHi - rcLo > bestSupport:
                                            bestSupport = rcHi - rcLo
                                            # bestSupport = hi - lo + rcHi - rcLo
                                            bestBase = base
                                    if bestSupport == 0:  #TODO: do this query more often (during ties?)
                                        for base in BASES:
                                            lo, hi = bwt.findIndicesOfStr(''.join(read[i+1:i+k]) + base)
                                            if hi - lo > bestSupport:
                                                bestSupport = hi - lo
                                                bestBase = base
                                    if bestBase != read[i+k]:
                                        changeMade = True
                                        corrected[i+k] = True
                                        read[i+k] = bestBase
                                    if corrected[i+k]:
                                        kmersEnd = min(len(read), i+2*k)
                                        # newCounts, _ = bwt.countStrandedSeqMatches(''.join(read[i+1:kmersEnd]), k)
                                        newCounts, _ = bwt.countStrandedSeqMatches(
                                            reverseComplement(''.join(read[i+1:kmersEnd])), k)
                                        newCounts = np.flipud(newCounts)
                                        trusted[i+1:kmersEnd-k+1] = newCounts > 0
                                        if False in trusted[i+1:kmersEnd-k+1]:
                                            trusted[i+1:kmersEnd-k+1] |= newCounts > hiThresh
                                            counts[i+1:kmersEnd-k+1] = newCounts
                            elif trusted[i+1] and not corrected[i]:  # err at read[i]
                                bestSupport = 0
                                newLo, newHi = bwt.findIndicesOfStr(''.join(read[i+1:i+k]))
                                for base in BASES:
                                    lo = bwt.getOccurrenceOfCharAtIndex(BASE_NUMS[base], newLo)
                                    hi = bwt.getOccurrenceOfCharAtIndex(BASE_NUMS[base], newLo)
                                    # lo, hi = bwt.findIndicesOfStr(base + ''.join(read[i+1:i+k]))
                                    # if hi - lo + rcHi - rcLo > bestSupport:
                                    if hi - lo > bestSupport:
                                        bestSupport = hi - lo
                                        # bestSupport = hi - lo + rcHi - rcLo
                                        bestBase = base
                                if bestSupport == 0:
                                    for base in BASES:
                                        rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(base + ''.join(read[i+1:i+k])))
                                        if rcHi - rcLo > bestSupport:
                                            bestSupport = rcHi - rcLo
                                            bestBase = base
                                if bestBase != read[i]:
                                    changeMade = True
                                    corrected[i] = True
                                    read[i] = bestBase
                                if corrected[i]:
                                    kmersStart = max(0, i-k)
                                    newCounts, _ = bwt.countStrandedSeqMatches(
                                        reverseComplement(''.join(read[kmersStart:i+k])), k)
                                    newCounts = np.flipud(newCounts)
                                    trusted[kmersStart:i+1] = newCounts > 0
                                    # counts[kmersStart:i+1] = newCounts
                                    if False in trusted[kmersStart:i+1]:
                                        newCounts, _ = bwt.countStrandedSeqMatches(''.join(read[kmersStart:i+k]), k)
                                        trusted[kmersStart:i+1] |= newCounts > hiThresh
            tmpFile.write(readName)
            tmpFile.write(''.join(read) + '\n')
            tmpFile.write(plus)
            tmpFile.write(qual)
    logging.info('Process %d finished in %.2f s', processNum, clock() - begin)
    return tmpFile.name


def willBeMain(inFilename, bwtDir, k=25, hiThresh=5, outFilename='corrected.fastq', numProcesses=1):
    logging.basicConfig(level=logging.DEBUG)
    begin = clock()
    pool = Pool(numProcesses)
    mapFunc = partial(correct, inFilename, bwtDir, k, hiThresh, numProcesses)
    fastqs = pool.map(mapFunc, range(numProcesses))
    logging.info('All child processes finished in %.2f s', clock() - begin)
    logging.info('Combining temporary files...')
    with open(outFilename, 'w') as outFile:
        for fastq in fastqs:
            with open(fastq) as inFile:
                for line in inFile:
                    outFile.write(line)
            remove(fastq)
    logging.info('Corrected reads saved in %s', outFilename)


def main():
    willBeMain('/playpen/sgreens/ecoli/EAS20_8/cov20.txt', '/playpen/sgreens/ecoli/msbwt20/rle_bwt/',
               outFilename='/playpen/sgreens/ecoli/msbwt20/corrected.fastq', numProcesses=4)
