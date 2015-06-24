from MUSCython import MultiStringBWTCython as msbwt, LCPGen
from MUSCython.BasicBWT cimport BasicBWT
from time import clock
from tempfile import NamedTemporaryFile
import logging
from multiprocessing import Pool
from os import remove, path
from string import maketrans
from functools import partial
import cython
import numpy as np
cimport numpy as np

cdef list BASES = ['A', 'C', 'G', 'T']

BASE_NUMS = {'A': 1, 'C': 2, 'G': 3, 'T': 5}
NUM_TO_BASE = ['A', 'C', 'G', 'N', 'T', '$']

TRAN_TAB = maketrans('ACGT', 'TGCA')

cdef bytes reverseComplement(seq):
    return seq[::-1].translate(TRAN_TAB)


def fastaParser(fasta, int start, int end):
    cdef long _
    with open(fasta, 'r') as fp:
        for _ in xrange(start*4):
                fp.next()
        for _ in xrange(start, end):
            yield fp.next(), fp.next(), fp.next(), fp.next()

@cython.boundscheck(False)
def correct(inFastqFile, bwtDir, k, hiThresh, numProcesses, processNum):
    begin = clock()
    cdef BasicBWT bwt = msbwt.loadBWT(bwtDir, False)
    cdef np.ndarray[np.uint8_t] corrected
    cdef np.ndarray[np.uint8_t] trusted
    cdef np.ndarray[np.uint64_t] counts
    cdef np.ndarray[np.uint64_t] revCounts
    cdef np.ndarray[np.uint64_t] newCounts
    cdef int readsPerProcess = (bwt.getSymbolCount(0) / numProcesses) + 1
    cdef int kmersEnd, kmersStart
    cdef long lo, hi, rcLo, rcHi, bestSupport, support
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
            read = list(origRead[:len(origRead)-1])
            revCounts = bwt.countStrandedSeqMatchesNoOther(reverseComplement(origRead[:len(origRead)-1]), k)
            revCounts = np.flipud(revCounts)
            trusted = (revCounts > 0).astype(np.uint8)
            corrected = np.zeros(len(read), dtype=np.uint8)
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
                                    newLo, newHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[i+1:i+k])))
                                    for base in BASES:
                                        rcLo = bwt.getOccurrenceOfCharAtIndex(BASE_NUMS[base.translate(TRAN_TAB)], newLo)
                                        rcHi = bwt.getOccurrenceOfCharAtIndex(BASE_NUMS[base.translate(TRAN_TAB)], newHi)
                                        lo, hi = bwt.findIndicesOfStr(<bytes> ''.join(read[i+1:i+k]) + base)
                                        support = hi - lo + rcHi - rcLo
                                        if support > bestSupport:
                                            bestSupport = support
                                            bestBase = base
                                    if bestBase != read[i+k]:
                                        changeMade = True
                                        corrected[i+k] = True
                                        read[i+k] = bestBase
                                        kmersEnd = min(len(read), i+2*k)
                                        newCounts = bwt.countStrandedSeqMatchesNoOther(
                                            reverseComplement(''.join(read[i+1:kmersEnd])), k)
                                        newCounts = np.flipud(newCounts)
                                        trusted[i+1:kmersEnd-k+1] = newCounts > 0
                                        if False in trusted[i+1:kmersEnd-k+1]:
                                            newCounts = bwt.countStrandedSeqMatchesNoOther(<bytes> ''.join(read[i+1:kmersEnd]), k)
                                            trusted[i+1:kmersEnd-k+1] |= newCounts > hiThresh
                            elif trusted[i+1] and not corrected[i]:  # err at read[i]
                                bestSupport = 0
                                newLo, newHi = bwt.findIndicesOfStr(<bytes> ''.join(read[i+1:i+k]))
                                for base in BASES:
                                    lo = bwt.getOccurrenceOfCharAtIndex(BASE_NUMS[base], newLo)
                                    hi = bwt.getOccurrenceOfCharAtIndex(BASE_NUMS[base], newHi)
                                    rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(base + ''.join(read[i+1:i+k])))
                                    support = hi - lo + rcHi - rcLo
                                    if support > bestSupport:
                                        bestSupport = support
                                        bestBase = base
                                if bestBase != read[i]:
                                    changeMade = True
                                    corrected[i] = True
                                    read[i] = bestBase
                                    kmersStart = max(0, i-k)
                                    newCounts = bwt.countStrandedSeqMatchesNoOther(
                                        reverseComplement(''.join(read[kmersStart:i+k])), k)
                                    newCounts = np.flipud(newCounts)
                                    trusted[kmersStart:i+1] = newCounts > 0
                                    if False in trusted[kmersStart:i+1]:
                                        newCounts = bwt.countStrandedSeqMatchesNoOther(<bytes> ''.join(read[kmersStart:i+k]), k)
                                        trusted[kmersStart:i+1] |= newCounts > hiThresh
            tmpFile.write(readName)
            tmpFile.write(''.join(read) + '\n')
            tmpFile.write(plus)
            tmpFile.write(qual)
    logging.info('Process %d finished in %.2f s', processNum, clock() - begin)
    return tmpFile.name


def willBeMain(inFilename, bwtDir, maxReadLen, k=25, hiThresh=5, outFilename='corrected.fastq', numProcesses=1):
    if not path.isfile(bwtDir+'lcps.npy'):
        buildLCP(bwtDir, maxReadLen)
    logging.basicConfig(level=logging.DEBUG)
    begin = clock()
    # correct(inFilename, bwtDir, k, hiThresh, 1, 0)
    # exit()
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


def buildLCP(bwtDir, maxReadLen):
    logging.basicConfig(level=logging.DEBUG)
    logging.info('Building LCP array')
    begin = clock()
    cdef np.ndarray[np.uint8_t] lcps = LCPGen.lcpGenerator(bwtDir, maxReadLen+1, logging.getLogger())
    np.save(bwtDir + 'lcps.npy', lcps)
    logging.info('Finished building LCP in %d s', clock() - begin)


def main():
    # buildLCP('/playpen/sgreens/fq_celegans/msbwt/bwt/', 101)
    willBeMain('/playpen/sgreens/ecoli/EAS20_8/cov20.txt', '/playpen/sgreens/ecoli/msbwt20/rle_bwt/', 101,
               outFilename='/playpen/sgreens/ecoli/msbwt20/corrected.fastq', numProcesses=4)
    # prefix = '/playpen/sgreens/fq_celegans/'
    # willBeMain(prefix + 'cov10.txt', prefix + '/msbwt/bwt/', 101,
    #            outFilename=prefix+'/msbwt/corrected.fastq', numProcesses=4)
