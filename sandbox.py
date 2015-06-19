from MUSCython import MultiStringBWTCython as msbwt, LCPGen
from logging import getLogger
from time import clock
from itertools import izip, izip_longest
import pysam
from collections import Counter
import cProfile
from random import randint
from string import maketrans
import numpy as np
import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plot
# cimport numpy as np

# cdef int K = 19
# cdef list BASES = ['A', 'C', 'G', 'T']
# cdef int READ_LEN = 100
# PATH = '/playpen/sgreens/ecoli/'
# BWT_PATH = '/playpen/sgreens/ecoli/bwt/'
# PATH = '/playpen/sgreens/fake/'
# BWT_PATH = '/playpen/sgreens/fake/bwt/'
PATH = '/playpen/sgreens/ecoli/msbwt20/'
BWT_PATH = '/playpen/sgreens/ecoli/msbwt20/rle_bwt2/'


K = 25
BASES = ['A', 'C', 'G', 'T']
BASE_TO_NUM = {'A': 1, 'C': 2, 'G': 3, 'T': 5}
NUM_TO_BASE = ['$', 'A', 'C', 'G', 'N', 'T']
READ_LEN = 100
THRESH = 1
HI_THRESH = 5
TRAN_TAB = maketrans('ACGT', 'TGCA')


def reverseComplement(seq):
    return seq[::-1].translate(TRAN_TAB)


def grouper(iterable, n, fillvalue=None):
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


# @profile
def correct():
    begin = clock()
    bwt = msbwt.loadBWT(BWT_PATH, True)
    trusted = np.zeros(READ_LEN - K + 1, dtype=np.uint8)
    print bwt.getSymbolCount(0), 'reads'
    with open(PATH + 'corrected2.fa', 'w') as fp:
        for readID in xrange(0, bwt.getSymbolCount(0)):
        # for readID in xrange(0, 1000):
            trusted.fill(0)
            origRead = bwt.recoverString(readID)
            read = list(origRead[1:])
            changeMade = True
            while changeMade:
                changeMade = False
                if read[-1] != 'N':
                    lo, hi = bwt.findIndicesOfStr(''.join(read[-K:]))
                    if hi - lo > HI_THRESH:
                        trusted[-1] = True
                    else:
                        rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[-K:])))
                        trusted[-1] = rcHi != rcLo and hi - lo + rcHi - rcLo > THRESH
                for start, end in izip(xrange(READ_LEN - K - 1, -1, -1), xrange(READ_LEN - 1, K - 1, -1)):
                    if not trusted[start]:
                        if read[start] != 'N':
                            query = ''.join(read[start:end])
                            lo, hi = bwt.findIndicesOfStr(query)
                            if hi - lo > HI_THRESH:
                                trusted[start] = True
                            else:
                                rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[start:end])))
                                trusted[start] = (rcHi != rcLo and hi - lo + rcHi - rcLo > THRESH)
                        if not trusted[start] and trusted[start + 1]:
                            # first base of k-mer is error
                            bestSupport = 0
                            lo, hi = bwt.findIndicesOfStr(''.join(read[start + 1:end]))
                            for newBase in BASES:
                                newLo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[newBase], lo)
                                newHi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[newBase], hi)
                                rcLo, rcHi = bwt.findIndicesOfStr(
                                    reverseComplement(newBase + ''.join(read[start + 1:end])))
                                if newHi - newLo + rcHi - rcLo > bestSupport:
                                    bestSupport = newHi - newLo + rcHi - rcLo
                                    if read[start] != newBase:
                                        changeMade = True
                                    read[start] = newBase
                                    trusted[start] = bestSupport > THRESH
                        elif trusted[start] and not trusted[start + 1]:
                            # base following this k-mer is error
                            bestSupport = 0
                            rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[start + 1:end])))
                            for newBase in BASES:
                                newLo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[reverseComplement(newBase)], rcLo)
                                newHi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[reverseComplement(newBase)], rcHi)
                                lo, hi = bwt.findIndicesOfStr(''.join(read[start + 1:end]) + newBase)
                                if newHi - newLo + hi - lo > bestSupport:
                                    bestSupport = newHi - newLo + hi - lo
                                    if read[end] != newBase:
                                        changeMade = True
                                    read[end] = newBase
            fp.write('>' + str(readID) + '\n')
            fp.write(''.join(read) + '\n')
    print 'Took', clock() - begin, 'seconds'


def lcpCorrect():
    begin = clock()
    bwt = msbwt.loadBWT(BWT_PATH, True)
    trusted = np.zeros(READ_LEN - K + 1, dtype=np.uint8)
    changed = np.empty(READ_LEN, dtype=np.uint8)
    trustedIndices = np.load(BWT_PATH + 'trustedIndices.npy')
    print bwt.getSymbolCount(0), 'reads'
    with open(PATH + 'corrected2.fa', 'w') as fp:
        # for readID in xrange(readsPerWorker * workerNum, min(readsPerWorker * (workerNum+1), bwt.getSymbolCount(0))):
        for readID in xrange(bwt.getSymbolCount(0)):
            trusted.fill(0)
            changed.fill(0)
            origRead, indices = bwt.recoverString(readID, True)
            read = list(origRead[1:])
            changeMade = True
            while changeMade:
                changeMade = False
                if read[-1] != 'N':  # untrusted if k-mer starts with an N
                    # if this k-mer is unchanged, see if high thresh is met
                    if True not in changed[READ_LEN - K:] and trustedIndices[indices[READ_LEN - K + 1]]:
                        trusted[-1] = True
                    else:  # k-mer is either changed or didn't meet high thresh
                        lo, hi = bwt.findIndicesOfStr(''.join(read[-K:]))
                        if hi - lo > HI_THRESH:
                            trusted[-1] = True
                        else:
                            rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[-K:])))
                            trusted[-1] = rcHi != rcLo and hi - lo + rcHi - rcLo > THRESH
                for start, end in izip(xrange(READ_LEN - K - 1, -1, -1), xrange(READ_LEN - 1, K - 1, -1)):
                    if not trusted[start]:
                        if read[start] != 'N':
                            if True not in changed[start:end] and trustedIndices[indices[start + 1]]:
                                trusted[start] = True
                            else:
                                lo, hi = bwt.findIndicesOfStr(''.join(read[start:end]))
                                if hi - lo > HI_THRESH:
                                    trusted[start] = True
                                else:
                                    rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[start:end])))
                                    trusted[start] = (rcHi != rcLo and hi - lo + rcHi - rcLo > THRESH)
                        if not trusted[start] and trusted[start + 1]:
                            # first base of k-mer is error
                            bestSupport = 0
                            lo, hi = bwt.findIndicesOfStr(''.join(read[start + 1:end]))
                            for newBase in BASES:
                                newLo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[newBase], lo)
                                newHi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[newBase], hi)
                                rcLo, rcHi = bwt.findIndicesOfStr(
                                    reverseComplement(newBase + ''.join(read[start + 1:end])))
                                if newHi - newLo + rcHi - rcLo > bestSupport:
                                    bestSupport = newHi - newLo + rcHi - rcLo
                                    if read[start] != newBase:
                                        changeMade = True
                                        changed[start] = True
                                    read[start] = newBase
                                    trusted[start] = bestSupport > THRESH
                        elif trusted[start] and not trusted[start + 1]:
                            # base following this k-mer is error
                            bestSupport = 0
                            rcLo, rcHi = bwt.findIndicesOfStr(reverseComplement(''.join(read[start + 1:end])))
                            for newBase in BASES:
                                newLo = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[reverseComplement(newBase)], rcLo)
                                newHi = bwt.getOccurrenceOfCharAtIndex(BASE_TO_NUM[reverseComplement(newBase)], rcHi)
                                lo, hi = bwt.findIndicesOfStr(''.join(read[start + 1:end]) + newBase)
                                if newHi - newLo + hi - lo > bestSupport:
                                    bestSupport = newHi - newLo + hi - lo
                                    if read[end] != newBase:
                                        changeMade = True
                                        changed[end] = True
                                    read[end] = newBase
            fp.write('>' + str(readID) + '\n')
            fp.write(''.join(read) + '\n')
    print 'Took', clock() - begin, 'seconds'


def superCorrect():
    begin = clock()
    bwt = msbwt.loadBWT(BWT_PATH, False)
    trusted = np.empty(READ_LEN-K+1, dtype=np.uint8)
    corrected = np.empty(READ_LEN, dtype=np.uint8)
    with open(PATH + 'corrected2.fa', 'w') as fp:
        for readID in xrange(bwt.getSymbolCount(0)):
            trusted.fill(0)
            corrString = ['0'] * READ_LEN
            origRead = bwt.recoverString(readID)
            origRead = origRead[1:] + '\n'
            origRead = 'CGGTCGCGCTATACTTTAGATGCCCAGGTCGCTGCATCATGGGTAATGAAGAATAAGGCTGGATAAAGCGACGTTGTGTACCGTCACTTTCTTCAATCGT\n'
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
                                        corrString[i+K] = '1'
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
                                    corrString[i] = '1'
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
            print ''.join(read)
            exit()
            fp.write(''.join(read) + '\n')
            fp.write(origRead)
            fp.write(''.join(corrString) + '\n')
    print 'Took', clock() - begin, 'seconds'


def runLengthCorrect():
    begin = clock()
    bwt = msbwt.loadBWT(BWT_PATH, True)
    changed = np.empty(READ_LEN, dtype=np.uint8)
    print bwt.getSymbolCount(0), 'reads'
    DOLLAR = 0
    N = 4
    with open(PATH + 'corrected2.fa', 'w') as fp:
        for readID in xrange(bwt.getSymbolCount(0)):
            origRead, indices = bwt.recoverString(readID, True)
            read = list(origRead[1:])
            changed.fill(0)
            doQuery = False
            for charPos, origIndex in reversed(zip(range(len(indices)-1), indices[2:] + [indices[0]])):
                try:
                    c = bwt.getCharAtIndex(origIndex)
                    if charPos < 50 and N != c != DOLLAR:
                        aboveIndex = origIndex - 1
                        belowIndex = origIndex + 1
                        # if doQuery:
                        #     i = bwt.findIndicesOfStr(''.join(read[charPos+1:]))[0]
                        #     if i == origIndex:
                        #         doQuery = False
                        #     else:
                        #         aboveIndex = i
                        #         belowIndex = i + 1
                            # print i, bwt.findIndicesOfStr(''.join(read[charPos+1:]))[0], bwt.findIndicesOfStr(origRead[charPos+2:])[0]
                        charAbove = bwt.getCharAtIndex(aboveIndex)
                        charBelow = bwt.getCharAtIndex(belowIndex)
                        if c != charAbove and c != charBelow:
                            if charAbove == N or charAbove == DOLLAR:
                                if N != charBelow != DOLLAR:
                                    read[charPos] = NUM_TO_BASE[charBelow]
                                    changed[charPos] = True
                                    doQuery = True
                            elif charBelow == N or charBelow == DOLLAR:
                                read[charPos] = NUM_TO_BASE[charAbove]
                                changed[charPos] = True
                                doQuery = True
                            else:
                                # nextIndex = bwt.getOccurrenceOfCharAtIndex(c, i)
                                nextIndex = indices[charPos + 1]
                                if abs(bwt.getOccurrenceOfCharAtIndex(charAbove, aboveIndex) - nextIndex) > \
                                        abs(bwt.getOccurrenceOfCharAtIndex(charBelow, belowIndex) - nextIndex):
                                    read[charPos] = NUM_TO_BASE[charAbove]
                                    changed[charPos] = True
                                    doQuery = True
                                else:
                                    read[charPos] = NUM_TO_BASE[charBelow]
                                    changed[charPos] = True
                                    doQuery = True
                except OverflowError:
                    pass  # happens when we try to do getCharAtIndex -1
            fp.write(''.join(read) + '\n')
            fp.write(origRead[1:] + '\n')
            fp.write(''.join(['1' if c else '0' for c in changed]) + '\n')
    print 'Took', clock() - begin, 'seconds'


def correctKmer(bwt, read):
    for i in xrange(0, K):
        for newBase in BASES:
            newLo, newHi = bwt.findIndicesOfStr(''.join(read[:i] + [newBase] + read[i + 1:]))
            if newHi - newLo > THRESH:
                return newBase, i
    return None, None


def bruteForceCorrect(path=PATH):
    begin = clock()
    bwt = msbwt.loadBWT(BWT_PATH, False)
    trustedIndices = np.load(BWT_PATH + 'trustedIndices.npy')
    print bwt.getSymbolCount(0), 'reads'
    with open(path + 'corrected2.fasta', 'w') as fp:
        for readID in xrange(bwt.getSymbolCount(0)):
            # for readID in xrange(1000):
            trusted = [0] * (READ_LEN - K + 1)
            changed = [False] * READ_LEN
            # origRead = bwt.recoverString(readID)[1:]
            origRead, indices = bwt.recoverString(readID, True)
            origRead = origRead[1:]
            # indices = indices[1:]
            read = list(origRead)
            corrString = ['0'] * READ_LEN
            isChanging = True
            while isChanging:
                isChanging = False
                for start, end in izip(xrange(READ_LEN - K, -1, -1), xrange(READ_LEN, K - 1, -1)):
                    if not trusted[start]:
                        if True not in changed[start:end]:
                            trusted[start] = trustedIndices[indices[start + 1]]
                        else:
                            lo, hi = bwt.findIndicesOfStr(''.join(read[start:end]))
                            trusted[start] = (hi - lo > THRESH)
                        if not trusted[start]:
                            newBase, pos = correctKmer(bwt, read[start:end])
                            if newBase:
                                pos += start
                                read[pos] = newBase
                                corrString[pos] = '1'
                                changed[pos] = True
                                trusted[start] = 1
                                isChanging = True
            fp.write(''.join(read) + '\n')
            fp.write(origRead + '\n')
            fp.write(''.join(corrString) + '\n')
    print 'Took', clock() - begin, 'seconds'
    print path


def findTrustedIndices():
    begin = clock()
    lcps = np.load(BWT_PATH + 'lcp.npy')
    # lcps = [randint(1, 60) for _ in xrange(100)]
    trustedIndices = np.empty(len(lcps), dtype=np.uint8)
    blockStart = 0
    for i, lcp in enumerate(lcps, 1):
        if lcp < K:
            trustedIndices[blockStart:i] = i - blockStart > THRESH
            blockStart = i
    # for i in xrange(100):
    # print lcps[i], trustedIndices[i]
    np.save(BWT_PATH + 'trustedIndices.npy', trustedIndices)
    print 'findTrustedIndices:', clock() - begin, 'seconds'


def samplePaired():
    coverage = 20
    numReads = 28428648
    refLength = 4641652
    randReads = [randint(0, numReads - 1) for _ in xrange((coverage * refLength) / (2 * READ_LEN))]
    randReads.sort()
    print len(randReads), 'reads'
    for fileSuffix in ('1', '2'):
        with open('/playpen/sgreens/ecoli/EAS20_8/paired_' + str(coverage) + fileSuffix + '.txt', 'w') as outFile, \
                open('/playpen/sgreens/ecoli/EAS20_8/s_6_' + fileSuffix + '.txt', 'r') as inFile:
            reads = inFile.read().split('\n')
            for i in randReads:
                outFile.writelines([line + '\n' for line in reads[i * 4:4 * (i + 1)]])


def sample():
    coverage = 30
    numReads = 28428648
    refLength = 4641652
    randReads = [randint(0, numReads - 1) for _ in xrange((coverage * refLength) / READ_LEN)]
    randReads.sort()
    split = np.searchsorted(randReads, numReads / 2)
    randReads1 = randReads[:split]
    randReads2 = randReads[split:]
    print len(randReads), 'reads'
    print len(randReads1), 'reads from first file'
    print len(randReads2), 'reads from second file'
    with open('/playpen/sgreens/ecoli/EAS20_8/cov' + str(coverage) + '.txt', 'w') as outFile:
        with open('/playpen/sgreens/ecoli/EAS20_8/s_6_1.txt', 'r') as inFile:
            reads = inFile.read().split('\n')
            for i in randReads1:
                outFile.writelines([line + '\n' for line in reads[i * 4:4 * (i + 1)]])
        with open('/playpen/sgreens/ecoli/EAS20_8/s_6_2.txt', 'r') as inFile:
            reads = inFile.read().split('\n')
            for i in [j - (numReads / 2) for j in randReads2]:
                outFile.writelines([line + '\n' for line in reads[i * 4:4 * (i + 1)]])


def countsOfSeq(bwt, seq, k=K):
    ca = []
    for start in xrange(0, READ_LEN - k + 1):
        lo, hi = bwt.findIndicesOfStr(seq[start:start + k])
        c = hi - lo
        if c >= 10:
            ca.append('+')
        else:
            ca.append(str(c))
    return ''.join(ca)


def pseudoRef(fileName):
    samfile = pysam.Samfile(fileName, 'rb')
    lineLength = 70
    with open('/playpen/sgreens/ecoli/sequence.fasta', 'r') as fp:
        rawRef = fp.read().split('\n')[1:]
        # ref = ''.join(fp.read().split('\n')[1:])
    for col in samfile.pileup():
        c = Counter()
        for read in col.pileups:
            if read.query_position:
                c[read.alignment.seq[read.query_position]] += 1
        lineNum = col.reference_pos / lineLength
        linePos = col.reference_pos % lineLength
        line = rawRef[lineNum]
        try:
            rawRef[lineNum] = line[:linePos] + c.most_common(1)[0][0] + line[linePos + 1:]
        except IndexError:
            pass
    with open('/playpen/sgreens/ecoli/pseudoRef.fasta', 'w') as fp:
        for line in rawRef:
            fp.write(line + '\n')


def summarizeBam(fileName):
    MATCH = 0
    DEL = 2
    INS = 1
    SOFT_CLIP = 4
    HARD_CLIP = 5
    samfile = pysam.Samfile(fileName, 'r')
    errs = 0
    bases = 0
    alignedBases = 0
    perfect = 0
    errTypes = [0, 0, 0]
    with open('/playpen/sgreens/ecoli/sequence.fasta', 'r') as fp:
        # with open('/playpen/sgreens/ecoli/pseudoRef.fasta', 'r') as fp:
        ref = ''.join(fp.read().split('\n')[1:])
    bwt = msbwt.loadBWT(BWT_PATH)
    for read in samfile.fetch():
        bases += READ_LEN
        refPos = read.pos
        readPos = 0
        mismatches = ['0'] * READ_LEN
        for op, count in read.cigar:
            if op != MATCH and op != SOFT_CLIP and op != HARD_CLIP:
                errs += 1
                errTypes[op] += 1
            elif op == MATCH:
                alignedBases += count
                if read.seq[readPos:readPos + count] != ref[refPos:refPos + count]:
                    for errPos, (a, b) in enumerate(zip(read.seq[readPos:readPos + count], ref[refPos:refPos + count])):
                        if a != b:
                            errs += 1
                            errTypes[MATCH] += 1
                            mismatches[readPos + errPos] = '1'
                elif len(read.cigar) == 1:
                    perfect += 1
            if op != DEL:
                readPos += count
            if op != INS and op != SOFT_CLIP:
                refPos += count
                if '1' in mismatches:
                    if read.is_reverse:
                        print read.qname, 'rc'
                        print reverseComplement(ref[read.pos:read.pos+READ_LEN])
                        print reverseComplement(read.seq)
                        print ''.join(mismatches)[::-1]
                        # print countsOfSeq(bwt, reverseComplement(read.seq))
                        # print countsOfSeq(bwt, read.seq)[::-1]
                        print countsOfSeq(bwt, reverseComplement(read.seq), 25)
                        print countsOfSeq(bwt, read.seq, 25)[::-1]
                    else:
                        print read.qname
                        print ref[read.pos:read.pos+READ_LEN]
                        print read.seq
                        print ''.join(mismatches)
                        # print countsOfSeq(bwt, read.seq)
                        # print countsOfSeq(bwt, reverseComplement(read.seq))[::-1]
                        print countsOfSeq(bwt, read.seq, 25)
                        print countsOfSeq(bwt, reverseComplement(read.seq), 25)[::-1]
                    print
    print bases, 'bases'
    print alignedBases, 'aligned bases'
    print errs, 'errors'
    print 'Error rate 1 per', bases / errs, 'bases'
    print perfect, 'perfect reads'
    print errTypes[0], 'Mismatches'
    print errTypes[1], 'Insertions'
    print errTypes[2], 'Deletions'


def compareQuals(fileName):
    MATCH = 0
    DEL = 2
    INS = 1
    SOFT_CLIP = 4
    HARD_CLIP = 5
    samfile = pysam.Samfile(fileName, 'r')
    bases = 0
    alignedBases = 0
    NUM_SAMPLES = 1000
    numErrs = 0
    errs = np.empty(NUM_SAMPLES, dtype=np.uint8)
    positions = np.empty(NUM_SAMPLES, dtype=np.uint8)
    quals = np.empty(NUM_SAMPLES, dtype=np.uint8)
    numSamples = 0
    with open('/playpen/sgreens/ecoli/sequence.fasta', 'r') as fp:
        ref = ''.join(fp.read().split('\n')[1:])
    for read in samfile.fetch():
        bases += READ_LEN
        refPos = read.pos
        readPos = 0
        for op, count in read.cigar:
            if op == MATCH:
                alignedBases += count
                if read.query_alignment_sequence[readPos:readPos + count] != ref[refPos:refPos + count]:
                    for errPos, (a, b) in enumerate(zip(read.seq[readPos:readPos + count], ref[refPos:refPos + count])):
                        if a != 'N':
                            if a != b:
                                errs[numSamples] = 1
                                positions[numSamples] = readPos + errPos
                                quals[numSamples] = read.query_qualities[readPos + errPos]
                                numSamples += 1
                                numErrs += 1
                            elif numErrs > numSamples / 2:  # if we need more non-error samples
                                errs[numSamples] = 0
                                positions[numSamples] = readPos + errPos
                                quals[numSamples] = read.query_qualities[readPos + errPos - read.query_alignment_start]
                                numSamples += 1
                            if numSamples >= NUM_SAMPLES:
                                plot.ylabel('Quality')
                                plot.xlabel('Position in read')
                                plot.scatter(positions[errs == 1], quals[errs == 1], c='r')
                                plot.scatter(positions[errs == 0], quals[errs == 0], c='b')
                                plot.title('Base quality vs position in read')
                                plot.figtext(.5, .85, 'Errors red, non-errors blue', ha='center')
                                plot.savefig('/csbiohome01/sgreens/Projects/refAlign/qual.png')
                                return
            if op != DEL:
                readPos += count
            if op != INS and op != SOFT_CLIP:
                refPos += count


def convertToFasta():
    readNum = 0
    with open('/playpen/sgreens/ecoli/corrected_combined_old.fa', 'r') as inFile:
        with open('/playpen/sgreens/ecoli/corrected_combined.fa', 'w') as outFile:
            for line in inFile:
                outFile.write('>' + str(readNum) + '\n')
                outFile.write(line)
                readNum += 1


def main():
    # np.save('/playpen/sgreens/fake/bwt/lcps.npy', LCPGen.lcpGenerator(BWT_PATH, READ_LEN+1, getLogger()))
    # cProfile.run('correct()')
    # compareQuals('/playpen/sgreens/ecoli/uncorrected20.sam')
    # correct()
    # runLengthCorrect()
    # superCorrect()
    # print 'msbwt';
    # summarizeBam('/playpen/sgreens/ecoli/msbwt20/msbwt.sam')
    summarizeBam('/playpen/sgreens/ecoli/msbwt20/msbwt_s.sam')
    # summarizeBam('/playpen/sgreens/ecoli/msbwtNoN20/msbwt.sam')
    # summarizeBam('/playpen/sgreens/ecoli/sga20/sga.sam')
    # print '\nsga';
    # summarizeBam('/playpen/sgreens/ecoli/msbwt20/rl.sam')
    # print '\nuncorrected';
    # summarizeBam('/playpen/sgreens/ecoli/uncorrected20.sam')
    # convertToFasta()


if __name__ == '__main__':
    main()