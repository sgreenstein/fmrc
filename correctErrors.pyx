#cython: boundscheck=False
#cython: wraparound=False

from MUSCython import MultiStringBWTCython as msbwt, LCPGen
from MUSCython.BasicBWT cimport BasicBWT
from time import clock
from tempfile import NamedTemporaryFile
import logging
import argparse
from multiprocessing import Pool
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
        lastPos = fp.tell()
        if fasta.endswith(('.fastq', '.fq')):
            while not fp.readline().startswith('@'):
                lastPos = fp.tell()
            fp.seek(lastPos)
            for _ in xrange(start*4):
                    fp.next()
            for _ in xrange(start, end):
                yield fp.next(), fp.next(), fp.next(), fp.next()
        elif fasta.endswith(('.fasta', '.fa')):
            while not fp.readline().startswith('>'):
                lastPos = fp.tell()
            fp.seek(lastPos)
            for _ in xrange(start*2):
                fp.next()
            for _ in xrange(start, end):
                yield fp.next(), fp.next(), '', ''
        else:
            raise Exception('%s does not have a correct fastq/fasta file extension' % fasta)


cpdef bytes correct(bytes inFile, bytes bwtDir, int k, int revCompThresh, int forwardThresh, int numProcesses, int processNum):
    begin = clock()
    cdef BasicBWT bwt = msbwt.loadBWT(bwtDir, False)
    cdef np.ndarray[np.uint8_t] corrected
    cdef np.ndarray[np.uint8_t] trusted
    cdef np.ndarray[np.uint64_t] counts
    cdef np.ndarray[np.uint64_t] revCounts
    cdef np.ndarray[np.uint64_t] newCounts
    cdef int readsPerProcess = (bwt.getSymbolCount(0) / numProcesses) + 1
    cdef int i, kmersEnd, kmersStart, readLen
    cdef long lo, hi, rcLo, rcHi, bestSupport, support
    cdef bytes readName, origRead, plus, qual, base
    cdef np.uint8_t bestBase
    cdef bint changeMade
    cdef char* read
    cdef bint printProgress = (processNum == numProcesses - 1)
    tmpFile = NamedTemporaryFile(suffix=inFile[inFile.rfind('.'):], delete=False)
    logging.info('Process %d / %d: correcting reads %d through %d', processNum, numProcesses,
                 readsPerProcess*processNum, min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))-1)
    cdef long readsDone = 0
    with tmpFile:
        for readName, origRead, plus, qual in \
                fastaParser(inFile, readsPerProcess*processNum,
                            min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))):
            read = origRead
            readLen = len(origRead)-1
            revCounts = bwt.countStrandedSeqMatchesNoOther((origRead[readLen-1::-1]).translate(TRAN_TAB), k)
            revCounts = np.flipud(revCounts)
            trusted = (revCounts > revCompThresh).astype(np.uint8)
            corrected = np.zeros(readLen, dtype=np.uint8)
            readsDone += 1
            if printProgress and not readsDone & 4194303:
                logging.info('%d reads done', readsPerProcess*processNum + readsDone)
            if False in trusted:
                counts = bwt.countStrandedSeqMatchesNoOther(origRead[:readLen], k)
                trusted |= counts > forwardThresh
                changeMade = True
                while changeMade:
                    changeMade = False
                    if False in trusted:
                        for i in xrange(len(trusted)-1):
                            if trusted[i] or read[i+k] == 'N':
                                if not trusted[i+1] or read[i+k] == 'N':  # err at read[i+k]
                                    bestSupport = 0
                                    bestBase = read[i+k]
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
                                        kmersEnd = min(readLen, i+2*k)
                                        newCounts = bwt.countStrandedSeqMatchesNoOther(
                                            read[kmersEnd-1:i:-1].translate(TRAN_TAB), k)
                                        newCounts = np.flipud(newCounts)
                                        trusted[i+1:kmersEnd-k+1] = newCounts > revCompThresh
                                        if False in trusted[i+1:kmersEnd-k+1]:
                                            newCounts = bwt.countStrandedSeqMatchesNoOther(read[i+1:kmersEnd], k)
                                            trusted[i+1:kmersEnd-k+1] |= newCounts > forwardThresh
                            elif trusted[i+1] and not corrected[i]:  # err at read[i]
                                bestSupport = 0
                                bestBase = read[i]
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
                                    trusted[kmersStart:i+1] = newCounts > revCompThresh
                                    if False in trusted[kmersStart:i+1]:
                                        newCounts = bwt.countStrandedSeqMatchesNoOther(read[kmersStart:i+k], k)
                                        trusted[kmersStart:i+1] |= newCounts > forwardThresh
            tmpFile.write(readName)
            tmpFile.write(read)
            tmpFile.write(plus)
            tmpFile.write(qual)
    logging.info('Process %d finished in %.2f s', processNum, clock() - begin)
    return tmpFile.name


def driver(inFilename, bwtDir, maxReadLen, k=25, revCompThresh=0, forwardThresh=5, outFilename='corrected.fastq', numProcesses=1):
    if not path.isfile(bwtDir+'lcps.npy'):
        buildLCP(bwtDir, maxReadLen)
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
    begin = clock()
    # correct(inFilename, bwtDir, k, revCompThresh, forwardThresh, 1, 0)
    # exit()
    pool = Pool(numProcesses)
    mapFunc = partial(correct, inFilename, bwtDir, k, revCompThresh, forwardThresh, numProcesses)
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
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
    logging.info('Building LCP array')
    begin = clock()
    cdef np.ndarray[np.uint8_t] lcps = LCPGen.lcpGenerator(bwtDir, maxReadLen+1, logging.getLogger())
    np.save(bwtDir + 'lcps.npy', lcps)
    logging.info('Finished building LCP in %d s', clock() - begin)


def main():
    # driver('/playpen/sgreens/ecoli/EAS20_8/cov20.fastq', '/playpen/sgreens/ecoli/msbwt20/rle_bwt/', 101,
    #            outFilename='/playpen/sgreens/ecoli/msbwt20/corrected.fastq', numProcesses=4)
    prefix = '/playpen/sgreens/fq_celegans/'
    driver(prefix + 'srr065388.fastq', prefix + '/msbwt60/bwt/', 101,
               outFilename=prefix+'/msbwt60/corrected.fastq', numProcesses=4)
    # parser = argparse.ArgumentParser()
    # parser.add_argument('inFile', help='the input fasta/fastq file of reads')
    # parser.add_argument('-b', '--bwtDir', metavar='bwtDir', dest='bwtDir', required=True, help='the directory containing the MSBWT')
    # parser.add_argument('-l', '--maxReadLength', metavar='maxReadLength', required=True, type=int, dest='maxReadLen',
    #                     help='the length of the longest read in the input file')
    # parser.add_argument('-k', default=25, metavar='k', type=int, dest='k', help='k-mer size')
    # parser.add_argument('-r', '--revCompThresh', metavar='revCompThresh', type=int, default=0, dest='revCompThresh',
    #                     help='k-mers whose reverse complement occurs more than this threshold will be trusted')
    # parser.add_argument('-f', '--forwardThresh', metavar='forwardThresh', type=int, default=5, dest='forwardThresh',
    #                     help='k-mers that occur more than this threshold will be trusted')
    # parser.add_argument('-o', metavar='outFile', dest='outFile', default='corrected',
    #                     help='the output fasta/fastq file of corrected reads')
    # parser.add_argument('-p', default=1, metavar='numProcesses', type=int, dest='numProcesses',
    #                     help='number of processes to use')
    # args = parser.parse_args()
    # fasta = ('.fa', '.fasta')
    # fastq = ('.fq', '.fastq')
    # outFile = args.outFile
    # if args.inFile.endswith(fasta):
    #     if outFile == 'corrected':
    #         outFile += '.fasta'
    # elif args.inFile.endswith(fastq):
    #     if outFile == 'corrected':
    #         outFile += '.fastq'
    # else:
    #     raise Exception('%s: Does not have correct fasta/fastq file extension' % args.inFile)
    # if not path.isfile(args.inFile):
    #     raise Exception('%s: No such file' % args.inFile)
    # if not 0 < args.k < args.maxReadLen:
    #     raise Exception('k must satisfy 0 < k < maxReadLen')
    # if args.forwardThresh < 1:
    #     raise Exception('forwardThresh must be a positive integer')
    # if args.revCompThresh < 0:
    #     raise Exception('revCompThresh must be a non-negative integer')
    # if args.maxReadLen <= 0:
    #     raise Exception('maxReadLength must be a positive integer')
    # try:
    #     msbwt.loadBWT(args.bwtDir)
    # except AttributeError:
    #     raise Exception('%s: Not a valid MSBWT directory' % args.bwtDir)
    # driver(args.inFile, args.bwtDir, args.maxReadLen, args.k, args.revCompThresh, args.forwardThresh, outFile,
    #        args.numProcesses)
