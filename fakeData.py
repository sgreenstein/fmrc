from random import randint, choice
from itertools import izip_longest
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plot

COVERAGE = 30
READ_LEN = 100
ERR_ODDS = 100  # i.e. 1 in ERR_ODDS chance of error
INDEL_ODDS = 1000
BASES = ['A', 'C', 'G', 'T']


def grouper(iterable, n, fillvalue=None):
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def createReads():
    with open('/playpen/sgreens/fake/ref.fasta', 'r') as fp:
        ref = ''.join(fp.read().split())  # remove newlines
    with open('/playpen/sgreens/fake/reads.txt', 'w') as fp:
        for readNum in xrange((COVERAGE * len(ref)) / READ_LEN):
            start = randint(0, len(ref) - READ_LEN)
            read = ref[start:start + READ_LEN]
            qualString = ''
            # add in errors
            # for pos in xrange(READ_LEN):
                # if randint(1, INDEL_ODDS) == 1:
                #
                # elif randint(1, ERR_ODDS) == 1:
                #     newBase = choice(BASES)
                #     while newBase == read[pos]:
                #         newBase = choice(BASES)
                #     read = read[:pos] + newBase + read[pos+1:]
                #     qualString += '1'  # error
                # else:
                #     qualString += '0'  # not error
            fp.write('@' + str(readNum) + '\n')
            fp.write(read + '\n')
            fp.write('+\n')
            fp.write(qualString + '\n')


def measurePerformance():
    with open('/playpen/sgreens/fake/ref.fasta', 'r') as fp:
        ref = ''.join(fp.read().split())  # remove newlines
    matches = 0
    with open('/playpen/sgreens/fake/reads.txt', 'r') as fp:
        for i, read in enumerate(fp):
            if (i-1) % 4 == 0:
                if read[:-1] in ref:
                    matches += 1
    print 'Before correction:', matches, 'reads matched', 'out of', i/4, 'total'
    matches = 0
    with open('/playpen/sgreens/fake/corrected2.fasta', 'r') as fp:
        for read, _, _ in grouper(fp, 3):
            # if i%1000 == 0:
            #     print i, 'reads done so far'
            if read[:-1] in ref:
                matches += 1
    print 'After correction:', matches, 'reads matched', 'out of', i/4, 'total'
    matches = 0
    with open('/playpen/sgreens/fake/reads.ec.fa', 'r') as fp:
        for _, read, _, _ in grouper(fp, 4):
            # if i%1000 == 0:
            #     print i, 'reads done so far'
            if read[:-1] in ref:
                matches += 1
    print 'SGA correction:', matches, 'reads matched', 'out of', i/4, 'total'


def correction():
    reads = {}
    with open('/playpen/sgreens/fake/reads.txt', 'r') as fp:
        for _, read, _, errString in grouper(fp, 4):
            reads[read[:-1]] = errString
    falsePos = 0
    falseNeg = 0
    truePos = 0
    trueNeg = 0
    missedErrPos = [0] * 100
    foundErrPos = [0] * 100
    errsInReadWhenMissed = [0] * 100
    errsInRead = [0] * 100
    with open('/playpen/sgreens/fake/corrected2.fa', 'r') as fp:
        for read, origRead, corrString in grouper(fp, 3):
            oldfp = falsePos
            errString = reads[origRead[:-1]][:-1]
            corrString = corrString[:-1]
            readErrs = errString.count('1')
            errsInRead[readErrs] += 1
            missedErr = False
            for pos, (err, corr) in enumerate(izip_longest(errString, corrString)):
                if err == '0':
                    if corr == '0':
                        trueNeg += 1
                    else:
                        falsePos += 1
                else:
                    if corr == '0':
                        falseNeg += 1
                        missedErrPos[pos] += 1
                        missedErr = True
                    else:
                        truePos += 1
                        foundErrPos[pos] += 1
            # if falsePos > oldfp:
            #     print origRead
            #     print 'err ', errString
            #     print 'corr', corrString
            #     print 'orig', origRead
            #     print 'new ', read
            #     exit()
            if missedErr:
                print origRead[:-1]
                print read[:-1]
                print errString
                print corrString
                print
                errsInReadWhenMissed[readErrs] += 1
    print 'Corrected %.3f%% of errors' % ((100.*truePos) / (truePos + falseNeg))
    print 'Corrected %.3f%% of non-errors' % ((100.*falsePos) / (trueNeg + falsePos))
    plot.bar(range(100), missedErrPos)
    plot.savefig('missed_err_pos.png')
    plot.clf()
    plot.bar(range(100), foundErrPos)
    plot.savefig('found_err_pos.png')
    plot.clf()
    plot.bar(range(100), errsInRead)
    plot.bar(range(100), errsInReadWhenMissed, color='r')
    plot.savefig('errs_in_read.png')
    print 'TP', truePos
    print 'TN', trueNeg
    print 'FP', falsePos
    print 'FN', falseNeg


def SGACorrection():
    falsePos = 0
    falseNeg = 0
    truePos = 0
    trueNeg = 0
    origReads = {}
    with open('/playpen/sgreens/fake/reads.txt', 'r') as fp:
        for readID, read, _, _ in grouper(fp, 4):
            readID = readID[:-1]
            read = read[:-1]
            origReads[readID] = read
    with open('/playpen/sgreens/fake/reads.ec.fa', 'r') as fp:
        for readID, read, _, errString in grouper(fp, 4):
            readID = readID[:-1]
            read = read[:-1]
            errString = errString[:-1]
            corrString = ['0' if a == b else '1' for a, b in izip_longest(origReads[readID], read)]
            for pos, (err, corr) in enumerate(izip_longest(errString, corrString)):
                if err == '0':
                    if corr == '0':
                        trueNeg += 1
                    else:
                        falsePos += 1
                else:
                    if corr == '0':
                        falseNeg += 1
                    else:
                        truePos += 1
    # print 'falsePos', falsePos
    # print 'falseNeg', falseNeg
    # print 'truePos', truePos
    # print 'trueNeg', trueNeg
    print 'Corrected %.3f%% of errors' % ((100.*truePos) / (truePos + falseNeg))
    print 'Corrected %.3f%% of non-errors' % ((100.*falsePos) / (trueNeg + falsePos))


def main():
    correction()
    # SGACorrection()
    # measurePerformance()


if __name__ == '__main__':
    main()