#cython: boundscheck=False
#cython: wraparound=False

from MUSCython import MultiStringBWTCython as msbwt, LCPGen
from MUSCython.BasicBWT cimport BasicBWT
from time import clock
from tempfile import NamedTemporaryFile
import logging
import argparse
from pysam import FastqFile
from multiprocessing import Pool
from os import remove, path
from string import maketrans
from functools import partial
import util
import cython
import numpy as np
cimport numpy as np

cdef list BASES = ['A', 'C', 'G', 'T']

# for converting into the msbwt's integer notation for bases
cdef dict BASE_TO_NUM = {'A': 1, 'C': 2, 'G': 3, 'T': 5}
cdef list NUM_TO_BASE = ['A', 'C', 'G', 'N', 'T', '$']

cdef bytes TRAN_TAB = maketrans('ACGT', 'TGCA')  # for reverse complementation

def fastaParser(bytes fasta, int start, int end):
    """ Generator for reading fasta/fastq files.
    Yields 4-tuple: seq id, seq, and: + and qual string if fastq, 2 empty strings if fasta
    :param fasta: path to the fasta or fastq file
    :param start: index of the read at which to start parsing
    :param end: index of the read at which to stop parsing (exclusive)
    :raise Exception: If input file contains neither a '@' nor a '>'
    """
    cdef long _
    cdef bytes fileType = None
    # determine file type by looking for sequence identifiers
    with open(fasta, 'r') as fp:
        lastPos = fp.tell()
        try:
            while not fileType:
                lastPos = fp.tell()
                line = fp.readline()
                if line.startswith('@'):
                    fileType = bytes('fastq')
                elif line.startswith('>'):
                    fileType = bytes('fasta')
        except StopIteration:
            raise Exception('Input file could not be parsed as fastq or fasta')
    # parse file
    fp = FastqFile(fasta)
    try:
        for _ in xrange(start):
            fp.next()
        if fileType == 'fastq':
            for _ in xrange(start, end):
                read = fp.next()
                yield '@' + read.name + '\n', read.sequence + '\n', '+\n', read.quality + '\n'
        else:
            for _ in xrange(start, end):
                read = fp.next()
                yield '>' + read.name + '\n', read.sequence + '\n', '', ''
    except StopIteration:
        logging.warning('Input %s file does not contain as many reads as the BWT', fileType)
    fp.close()


cpdef bytes correct(bytes inFile, bytes bwtDir, int k, int thresh, bint filterReads, int numProcesses, int processNum):
    """ Gets reads from a fasta/fastq file, corrects them, and writes a new fasta/fastq file
    :param inFile: path to the input fasta/fastq file
    :param bwtDir: path to the directory in which the msbwt is
    :param k: k-mer length
    :param thresh: if a k-mer's count exceeds this threshold, it is trusted
    :param numProcesses: how many concurrent processes there are
    :param processNum: the 0-based index of this process
    :return: path to the fasta/fastq file of corrected reads
    """
    begin = clock()
    cdef BasicBWT bwt = msbwt.loadBWT(bwtDir, True)
    cdef np.ndarray[np.uint8_t] corrected  # true if the base at that index has been corrected
    cdef np.ndarray[np.uint8_t] trusted  # stores trust of k-mers
    cdef np.ndarray[np.uint64_t] counts  # counts of forward k-mers
    cdef np.ndarray[np.uint64_t] revCounts  # counts of revcomp k-mers
    cdef np.ndarray[np.uint64_t] newCounts
    cdef int readsPerProcess = (bwt.getSymbolCount(0) / numProcesses) + 1
    cdef int i, kmersEnd, kmersStart, readLen
    cdef long lo, hi, rcLo, rcHi, bestSupport, support
    cdef bytes readName, origRead, plus, qual, base
    cdef np.uint8_t bestBase
    cdef bint changeMade
    cdef char* read
    cdef bint printProgress = (processNum == numProcesses - 1)
    cdef long readsDone = 0

    tmpFile = NamedTemporaryFile(dir=bwtDir, suffix=inFile[inFile.rfind('.'):], delete=False)
    logging.info('Process %d / %d: correcting reads %d through %d', processNum, numProcesses,
                 readsPerProcess*processNum, min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))-1)
    with tmpFile:
        for readName, origRead, plus, qual in \
                fastaParser(inFile, readsPerProcess*processNum,
                            min(readsPerProcess * (processNum+1), bwt.getSymbolCount(0))):
            # get counts for revcomp k-mers
            read = origRead
            readLen = len(origRead)-1
            revCounts = bwt.countStrandedSeqMatchesNoOther((origRead[readLen-1::-1]).translate(TRAN_TAB), k)
            revCounts = np.flipud(revCounts)
            trusted = (revCounts > thresh).astype(np.uint8)
            corrected = np.zeros(readLen, dtype=np.uint8)
            readsDone += 1
            if printProgress and not readsDone & 1048575:
                logging.info('%d reads done', readsPerProcess*processNum + readsDone)
            if False in trusted:
                # not all trusted from revcomp counts. Get forward counts
                counts = bwt.countStrandedSeqMatchesNoOther(origRead[:readLen], k)
                trusted |= counts > thresh + 1
                changeMade = True
                while changeMade:
                    # some k-mer(s) still untrusted based on forward and revcomp counts
                    # keep doing passes until we can't correct any more
                    changeMade = False
                    if False in trusted:
                        for i in xrange(len(trusted)-1):
                            if trusted[i] or read[i+k] == 'N':
                                if not trusted[i+1] or read[i+k] == 'N':
                                    # err at read[i+k]
                                    bestSupport = 0
                                    bestBase = read[i+k]
                                    newLo, newHi = bwt.findIndicesOfStr(read[i+k-1:i:-1].translate(TRAN_TAB))
                                    # try alternate bases
                                    for base in BASES:
                                        rcLo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[base.translate(TRAN_TAB)], newLo)
                                        rcHi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[base.translate(TRAN_TAB)], newHi)
                                        lo, hi = bwt.findIndicesOfStr(read[i+1:i+k] + base)
                                        support = hi - lo + rcHi - rcLo
                                        if support > bestSupport:
                                            bestSupport = support
                                            bestBase = ord(base)
                                    # if better base, get new counts and compute new trust
                                    if bestBase != read[i+k]:
                                        changeMade = True
                                        corrected[i+k] = True
                                        read[i+k] = bestBase
                                        kmersEnd = min(readLen, i+2*k)
                                        newCounts = bwt.countStrandedSeqMatchesNoOther(
                                            read[kmersEnd-1:i:-1].translate(TRAN_TAB), k)
                                        newCounts = np.flipud(newCounts)
                                        trusted[i+1:kmersEnd-k+1] = newCounts > thresh
                                        if False in trusted[i+1:kmersEnd-k+1]:
                                            newCounts = bwt.countStrandedSeqMatchesNoOther(read[i+1:kmersEnd], k)
                                            trusted[i+1:kmersEnd-k+1] |= newCounts > thresh
                            elif trusted[i+1] and not corrected[i]:
                                # err at read[i]
                                bestSupport = 0
                                bestBase = read[i]
                                newLo, newHi = bwt.findIndicesOfStr(read[i+1:i+k])
                                # try alternate bases
                                for base in BASES:
                                    lo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[base], newLo)
                                    hi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[base], newHi)
                                    rcLo, rcHi = bwt.findIndicesOfStr((read[i+k-1:i:-1] + base).translate(TRAN_TAB))
                                    support = hi - lo + rcHi - rcLo
                                    if support > bestSupport:
                                        bestSupport = support
                                        bestBase = ord(base)
                                # if better base, get new counts and compute new trust
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
                                    trusted[kmersStart:i+1] = newCounts > thresh
                                    if False in trusted[kmersStart:i+1]:
                                        newCounts = bwt.countStrandedSeqMatchesNoOther(read[kmersStart:i+k], k)
                                        trusted[kmersStart:i+1] |= newCounts > thresh
            if not filterReads or True in trusted:
                tmpFile.write(readName)
                tmpFile.write(read)
                tmpFile.write(plus)
                tmpFile.write(qual)
    logging.info('Process %d finished in %.2f s', processNum, clock() - begin)
    return tmpFile.name


def driver(inFilename, bwtDir, maxReadLen, k, thresh, filterReads, outFilename, numProcesses):
    """ Spawns processes that correct reads
    :param inFilename: path to the input fasta/fastq file containing reads
    :param bwtDir: path to the directory with the bwt in it
    :param maxReadLen: length of the longest read in the bwt
    :param k: k-mer length
    :param thresh: if a k-mer's count exceeds this threshold, it is trusted
    :param filterReads: if True, reads that have no trusted k-mers are not output
    :param outFilename: name of the fasta/fastq file to which to write corrected reads
    :param numProcesses: number of concurrent processes to use
    """
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    if not path.isfile(path.join(bwtDir, 'lcps.npy')):
        buildLCP(bwtDir, maxReadLen)
    begin = clock()
    pool = Pool(numProcesses)
    mapFunc = partial(correct, inFilename, bwtDir, k, thresh, filterReads, numProcesses)
    fastqs = pool.map(mapFunc, range(numProcesses))
    logging.info('All child processes finished in %.2f s', clock() - begin)
    begin = clock()
    logging.info('Combining temporary files...')
    with open(outFilename, 'w') as outFile:
        # copy headers from input file
        with open(inFilename) as inFile:
            read = inFile.readline()
            while not read.startswith(('@', '>')):
                outFile.write(read)
                read = inFile.readline()
        for fastq in fastqs:
            with open(fastq) as inFile:
                for line in inFile:
                    outFile.write(line)
            remove(fastq)
    logging.info('Combined temporary files in %.2f s', clock() - begin)
    logging.info('Corrected reads saved in %s', outFilename)


def buildLCP(bwtDir, maxReadLen):
    logging.info('Building LCP array')
    begin = clock()
    cdef np.ndarray[np.uint8_t] lcps = LCPGen.lcpGenerator(bwtDir, maxReadLen+1, logging.getLogger())
    np.save(bwtDir + 'lcps.npy', lcps)
    logging.info('Finished building LCP in %d s', clock() - begin)


def main():
    # set up parser
    parser = argparse.ArgumentParser()
    parser.add_argument('bwtDir', help='the directory containing the MSBWT')
    parser.add_argument('maxReadLength', type=int, help='the length of the longest read in the input file')
    parser.add_argument('inFile', help='the input fasta/fastq file of reads')
    parser.add_argument('-f', '--filter', dest='filterReads', action='store_true', default=False,
                        help='filter reads with uncorrectable errors')
    parser.add_argument('-k', default=25, metavar='k', type=int, dest='k', help='k-mer size')
    parser.add_argument('-t', '--threshold', metavar='threshold', type=int, default=2, dest='thresh',
                        help='k-mers that occur more than this threshold will be trusted')
    parser.add_argument('-o', metavar='outFile', dest='outFile', default='corrected',
                        help='the output fasta/fastq file of corrected reads')
    parser.add_argument('-p', default=1, metavar='numProcesses', type=int, dest='numProcesses',
                        help='number of processes to use')
    parser.add_argument('-v', '--version', action='version', version='fmrc version ' + util.VERSION)

    # parse and check args
    args = parser.parse_args()
    fastq = ('.fq', '.fastq')
    outFile = args.outFile
    if outFile == 'corrected':
        outFile += args.inFile[args.inFile.rfind('.'):]
    if not path.isfile(args.inFile):
        raise Exception('%s: No such file' % args.inFile)
    if not 0 < args.k < args.maxReadLength:
        raise Exception('k must satisfy 0 < k < maxReadLength')
    if args.thresh < 1:
        raise Exception('threshold must be a positive integer')
    if args.maxReadLength <= 0:
        raise Exception('maxReadLength must be a positive integer')
    try:
        msbwt.loadBWT(args.bwtDir)
    except AttributeError:
        raise Exception('%s: Not a valid MSBWT directory' % args.bwtDir)

    # do correction
    driver(args.inFile, args.bwtDir, args.maxReadLength, args.k, args.thresh, args.filterReads,
           outFile, args.numProcesses)
