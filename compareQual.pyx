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
PATH = '/playpen/sgreens/ecoli/msbwt20/'
BWT_PATH = '/playpen/sgreens/ecoli/msbwt20/rle_bwt/'

BASE_TO_NUM = {'A': 1, 'C': 2, 'G': 3, 'T': 5}
NUM_TO_BASE = ['A', 'C', 'G', 'N', 'T', '$']

TRAN_TAB = maketrans('ACGT', 'TGCA')

def reverseComplement(seq):
    return seq[::-1].translate(TRAN_TAB)


def grouper(iterable, n, fillvalue=None):
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def fastaParser(fasta, int start, int end):
    cdef long _
    with open(fasta, 'r') as fp:
        for _ in xrange(start*4):
                fp.next()
        for _ in xrange(start, end):
            yield fp.next(), fp.next(), fp.next(), fp.next()


def correct(wn_nw):
    begin = clock()
    workerNum = wn_nw[0]
    numWorkers = wn_nw[1]
    bwt = msbwt.loadBWT(BWT_PATH, True)
    cdef unsigned long lo, hi, newHi, newLo, rcLo, rcHi
    cdef int start, end, numTrusted, numUntrusted
    cdef np.ndarray[np.uint8_t] trusted = np.zeros(READ_LEN-K+1, dtype=np.uint8)
    cdef int readsPerWorker = (bwt.getSymbolCount(0) / (2*numWorkers)) + 1  # reads per worker per file
    print 'Worker' , workerNum, '/', numWorkers, 'doing', readsPerWorker*workerNum, 'through', min(readsPerWorker * (workerNum+1), bwt.getSymbolCount(0)/2)
    # for fileSuffix in ('1', '2'):
    fileSuffix = ''

    numSamples = 10000
    errs = np.empty(numSamples, dtype=np.uint8)
    numErrs = 0
    nonErrs = np.empty(numSamples, dtype=np.uint8)
    numNonErrs = 0

    with open(PATH + 'corrected' + fileSuffix + '_' + str(workerNum) + '.fastq', 'w') as fp:
        # for readID in xrange(readsPerWorker * workerNum, min(readsPerWorker * (workerNum+1), bwt.getSymbolCount(0))):
        for readName, origRead, plus, qual in\
            fastaParser('/playpen/sgreens/ecoli/EAS20_8/cov20' + fileSuffix + '.txt',
                        readsPerWorker*workerNum,
                        min(readsPerWorker * (workerNum+1), bwt.getSymbolCount(0)/2)):
            trusted.fill(0)
            read = list(origRead[:-1])
            changeMade = True

            corrections = np.zeros(READ_LEN, dtype=np.uint8)

            while changeMade:
                changeMade = False
                if read[-1] != 'N':
                    lo, hi = bwt.findIndicesOfStr(''.join(read[-K:]))
                    if hi - lo > HI_THRESH:
                        trusted[-1] = True
                    else:
                        rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[-K:])))
                        trusted[-1] = rcHi != rcLo and hi - lo + rcHi - rcLo > THRESH
                for start, end in izip(xrange(READ_LEN-K-1, -1, -1), xrange(READ_LEN-1, K-1, -1)):
                    if not trusted[start]:
                        if read[start] != 'N':
                            lo, hi = bwt.findIndicesOfStr(''.join(read[start:end]))
                            if hi - lo > HI_THRESH:
                                trusted[start] = True
                            else:
                                rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[start:end])))
                                trusted[start] = (rcHi != rcLo and hi - lo + rcHi - rcLo > THRESH)
                        if not trusted[start] and trusted[start+1]:
                            # first base of k-mer is error
                            bestSupport = 0
                            lo, hi = bwt.findIndicesOfStr(''.join(read[start+1:end]))
                            for newBase in BASES:
                                newLo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[newBase], lo)
                                newHi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[newBase], hi)
                                rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(newBase + ''.join(read[start+1:end])))
                                if newHi - newLo + rcHi - rcLo > bestSupport:
                                    bestSupport = newHi - newLo + rcHi - rcLo
                                    if read[start] != newBase:
                                        changeMade = True
                                        corrections[start] = 1
                                    read[start] = newBase
                                    trusted[start] = bestSupport > THRESH
                            if bestSupport <= THRESH:
                                corrections[start] = 2
                        elif trusted[start] and not trusted[start+1]:
                            # base following this k-mer is error
                            bestSupport = 0
                            rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[start+1:end])))
                            for newBase in BASES:
                                newLo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[reverseComplement(newBase)], rcLo)
                                newHi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[reverseComplement(newBase)], rcHi)
                                lo, hi = bwt.findIndicesOfStr(''.join(read[start+1:end]) + newBase)
                                if newHi - newLo + hi - lo > bestSupport:
                                    bestSupport = newHi - newLo + hi - lo
                                    if read[end] != newBase:
                                        changeMade = True
                                        # corrections[end] = 1
                                    read[end] = newBase
            if 1 in corrections:
                errs[numErrs] = ord(qual[choice(np.where(corrections == 1)[0])]) - ord('!')
                numErrs += 1
                nonErrs[numNonErrs] = ord(qual[choice(np.where(corrections == 0)[0])]) - ord('!')
                numNonErrs += 1
            # write corrected read to fastq
            fp.write(readName)
            fp.write(''.join(read) + '\n')
            fp.write(plus)
            fp.write(qual)
            if numErrs == numSamples:
                break
    plot.title('Distribution of base qualities')
    plot.xlabel('Base quality score')
    plot.ylabel('Number of bases')
    plot.hist(errs, bins=range(ord('~') - ord('!')), alpha=0.5, label='Error')
    plot.hist(nonErrs, bins=range(ord('~') - ord('!')), alpha=0.5, label='Non-error')
    plot.xlim(0, 45)
    plot.legend(loc='upper center')
    plot.savefig('/csbiohome01/sgreens/Projects/refAlign/qual_distrib.png')
    plot.ylim(0, 1000)
    plot.savefig('/csbiohome01/sgreens/Projects/refAlign/qual_low_distrib.png')
    print 'Took', clock() - begin, 'seconds'



def superCorrect(inFastqFile, bwtDir, int k, int loThresh, int hiThresh, int numProcesses, int processNum):
    begin = clock()
    bwt = msbwt.loadBWT(bwtDir, False)
    cdef np.ndarray[np.uint8_t] corrected
    cdef np.ndarray[np.uint8_t] trusted
    cdef np.ndarray[np.uint64_t] counts
    cdef np.ndarray[np.uint64_t] revCounts
    cdef np.ndarray[np.uint64_t] newCounts
    cdef int readsPerProcess = (bwt.getSymbolCount(0) / numProcesses) + 1
    cdef int i, kmersEnd, kmersStart
    cdef long lo, hi, rcLo, rcHi, bestSupport, support
    cdef bytes readName, origRead, plus, qual, base
    cdef np.uint8_t bestBase
    cdef bint changeMade
    cdef char* read
    numSamples = 10000
    errs = np.empty(numSamples, dtype=np.uint8)
    numErrs = 0
    nonErrs = np.empty(numSamples, dtype=np.uint8)
    numNonErrs = 0
    logging.info('Process %d / %d: correcting reads %d through %d', processNum, numProcesses,
                 readsPerProcess*processNum, min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))-1)
    cdef long readsDone = 0
    for readName, origRead, plus, qual in \
            fastaParser(inFastqFile, readsPerProcess*processNum,
                        min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))):
        read = origRead
        revCounts = bwt.countStrandedSeqMatchesNoOther((origRead[len(origRead)-2::-1]).translate(TRAN_TAB), k)
        revCounts = np.flipud(revCounts)
        trusted = (revCounts > loThresh).astype(np.uint8)
        corrected = np.zeros(len(origRead)-1, dtype=np.uint8)
        readsDone += 1
        # if not readsDone & 1023:
        #     logging.debug('Finished %d reads', readsDone)
        if False in trusted:
            counts = bwt.countStrandedSeqMatchesNoOther(origRead[:len(origRead)-1], k)
            trusted |= counts > hiThresh
            changeMade = True
            while changeMade:
                changeMade = False
                if False in trusted:
                    for i in xrange(len(trusted)-1):
                        if trusted[i]:
                            if numNonErrs < len(nonErrs):
                                nonErrs[numNonErrs] = ord(qual[i]) - ord('!')
                                numNonErrs += 1
                        if trusted[i] or read[i+k] == 'N':
                            if not trusted[i+1] or read[i+k] == 'N':  # err at read[i+k]
                                bestSupport = 0
                                newLo, newHi = bwt.findIndicesOfStr(read[i+k-1:i:-1].translate(TRAN_TAB))
                                for base in BASES:
                                    rcLo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[base.translate(TRAN_TAB)], newLo)
                                    rcHi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[base.translate(TRAN_TAB)], newHi)
                                    lo, hi = bwt.findIndicesOfStr(read[i+1:i+k] + base)
                                    support = hi - lo + rcHi - rcLo
                                    if support > bestSupport:
                                        bestSupport = support
                                        bestBase = ord(base)
                                if bestBase != read[i+k]:
                                    changeMade = True
                                    corrected[i+k] = True
                                    read[i+k] = bestBase
                                    kmersEnd = min(len(origRead)-1, i+2*k)
                                    newCounts = bwt.countStrandedSeqMatchesNoOther(
                                        read[kmersEnd-1:i:-1].translate(TRAN_TAB), k)
                                    newCounts = np.flipud(newCounts)
                                    trusted[i+1:kmersEnd-k+1] = newCounts > loThresh
                                    try:
                                        errs[numErrs] = ord(qual[i+k]) - ord('!')
                                        numErrs += 1
                                    except IndexError:
                                        return errs, nonErrs
                                    if False in trusted[i+1:kmersEnd-k+1]:
                                        newCounts = bwt.countStrandedSeqMatchesNoOther(read[i+1:kmersEnd], k)
                                        trusted[i+1:kmersEnd-k+1] |= newCounts > hiThresh
                        elif trusted[i+1] and not corrected[i]:  # err at read[i]
                            bestSupport = 0
                            newLo, newHi = bwt.findIndicesOfStr(read[i+1:i+k])
                            for base in BASES:
                                lo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[base], newLo)
                                hi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[base], newHi)
                                rcLo, rcHi = bwt.findIndicesOfStr((read[i+k-1:i:-1] + base).translate(TRAN_TAB))
                                support = hi - lo + rcHi - rcLo
                                if support > bestSupport:
                                    bestSupport = support
                                    bestBase = ord(base)
                            if bestBase != read[i]:
                                changeMade = True
                                corrected[i] = True
                                read[i] = bestBase
                                try:
                                    errs[numErrs] = ord(qual[i]) - ord('!')
                                    numErrs += 1
                                except IndexError:
                                    return errs, nonErrs
                                if i-k <= 0:
                                    kmersStart = 0
                                    newCounts = bwt.countStrandedSeqMatchesNoOther(
                                        read[i+k-1::-1].translate(TRAN_TAB), k)
                                else:
                                    kmersStart = i-k
                                    newCounts = bwt.countStrandedSeqMatchesNoOther(
                                        read[i+k-1:kmersStart-1:-1].translate(TRAN_TAB), k)
                                newCounts = np.flipud(newCounts)
                                trusted[kmersStart:i+1] = newCounts > loThresh
                                if False in trusted[kmersStart:i+1]:
                                    newCounts = bwt.countStrandedSeqMatchesNoOther(read[kmersStart:i+k], k)
                                    trusted[kmersStart:i+1] |= newCounts > hiThresh


def plotDistrib():
    errs, nonErrs = superCorrect('/playpen/sgreens/ecoli/EAS20_8/cov20.txt', '/playpen/sgreens/ecoli/msbwt20/bwt',
                                 25, 0, 5, 1, 0)
    plot.title('Distribution of base qualities')
    plot.xlabel('Base quality score')
    plot.ylabel('Number of bases')
    plot.hist(errs, bins=range(ord('~') - ord('!')), alpha=0.5, label='Error')
    plot.hist(nonErrs, bins=range(ord('~') - ord('!')), alpha=0.5, label='Non-error')
    plot.xlim(0, 45)
    plot.legend(loc='upper center')
    plot.savefig('/csbiohome01/sgreens/Projects/refAlign/qual_distrib.png')
    plot.ylim(0, 1000)
    plot.savefig('/csbiohome01/sgreens/Projects/refAlign/qual_low_distrib.png')

def main():
    correct((0, 1))
    exit()


if __name__ == '__main__':
    main()