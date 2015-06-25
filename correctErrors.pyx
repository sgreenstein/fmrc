#cython: boundscheck=False
#cython: wraparound=False

from MUSCython import MultiStringBWTCython as msbwt, LCPGen
from MUSCython.BasicBWT cimport BasicBWT
from time import clock
from tempfile import NamedTemporaryFile
import logging
from multiprocessing import Pool
from exceptions import OSError
from os import remove, path
from string import maketrans
from functools import partial
import cython
import numpy as np
cimport numpy as np

cdef list BASES = ['A', 'C', 'G', 'T']

cdef dict BASE_TO_NUM = {'A': 1, 'C': 2, 'G': 3, 'T': 5}
cdef list NUM_TO_BASE = ['A', 'C', 'G', 'N', 'T', '$']

cdef bytes TRAN_TAB = maketrans('ACGT', 'TGCA')

def fastaParser(bytes fasta, int start, int end):
    cdef long _
    with open(fasta, 'r') as fp:
        for _ in xrange(start*4):
                fp.next()
        for _ in xrange(start, end):
            yield fp.next(), fp.next(), fp.next(), fp.next()


cpdef bytes correct(bytes inFastqFile, bytes bwtDir, int k, int loThresh, int hiThresh, int numProcesses, int processNum):
    begin = clock()
    cdef BasicBWT bwt = msbwt.loadBWT(bwtDir, False)
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
    tmpFile = NamedTemporaryFile(suffix='.fastq', delete=False)
    logging.info('Process %d / %d: correcting reads %d through %d', processNum, numProcesses,
                 readsPerProcess*processNum, min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))-1)
    cdef long readsDone = 0
    with tmpFile:
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
            tmpFile.write(readName)
            tmpFile.write(read)
            tmpFile.write(plus)
            tmpFile.write(qual)
    logging.info('Process %d finished in %.2f s', processNum, clock() - begin)
    return tmpFile.name


def willBeMain(inFilename, bwtDir, maxReadLen, k=25, loThresh=0, hiThresh=5, outFilename='corrected.fastq', numProcesses=1):
    if not path.isfile(bwtDir+'lcps.npy'):
        buildLCP(bwtDir, maxReadLen)
    logging.basicConfig(level=logging.DEBUG)
    begin = clock()
    # correct(inFilename, bwtDir, k, loThresh, hiThresh, 1, 0)
    # exit()
    pool = Pool(numProcesses)
    mapFunc = partial(correct, inFilename, bwtDir, k, loThresh, hiThresh, numProcesses)
    fastqs = pool.map(mapFunc, range(numProcesses))
    logging.info('All child processes finished in %.2f s', clock() - begin)
    begin = clock()
    logging.info('Combining temporary files...')
    with open(outFilename, 'w') as outFile:
        for fastq in fastqs:
            with open(fastq) as inFile:
                for line in inFile:
                    outFile.write(line)
            remove(fastq)
    logging.info('Combined temporary files in %.2f s', clock() - begin)
    logging.info('Corrected reads saved in %s', outFilename)


def buildLCP(bwtDir, maxReadLen):
    logging.basicConfig(level=logging.DEBUG)
    logging.info('Building LCP array')
    begin = clock()
    cdef np.ndarray[np.uint8_t] lcps = LCPGen.lcpGenerator(bwtDir, maxReadLen+1, logging.getLogger())
    np.save(bwtDir + 'lcps.npy', lcps)
    logging.info('Finished building LCP in %d s', clock() - begin)


def main():
    # buildLCP('/playpen/sgreens/fq_celegans/msbwt/bwt/', 101)
    # willBeMain('/playpen/sgreens/ecoli/EAS20_8/cov20.txt', '/playpen/sgreens/ecoli/msbwt20/rle_bwt/', 101,
    #            outFilename='/playpen/sgreens/ecoli/msbwt20/corrected.fastq', numProcesses=4)
    prefix = '/playpen/sgreens/fq_celegans/'
    willBeMain(prefix + 'srr065388.fastq', prefix + '/msbwt60/bwt/', 101,
               outFilename=prefix+'/msbwt60/corrected.fastq', numProcesses=4)
