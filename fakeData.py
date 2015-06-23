from random import randint, choice
from itertools import izip_longest
from string import maketrans
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plot

COVERAGE = 20
READ_LEN = 100
ERR_ODDS = 100  # i.e. 1 in ERR_ODDS chance of error
INDEL_ODDS = 1000
BASES = ['A', 'C', 'G', 'T']

TRAN_TAB = maketrans('ACGT', 'TGCA')


def reverseComplement(seq):
    return seq[::-1].translate(TRAN_TAB)


def grouper(iterable, n, fillvalue=None):
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def createReads():
    with open('/playpen/sgreens/ecoli/sequence.fasta', 'r') as fp:
        ref = ''.join(fp.read().split('\n')[1:])  # remove newlines
    print ref[:100]
    ref = ref[:100000]
    with open('/playpen/sgreens/ecoli/EAS20_8/fake20.txt', 'w') as fp:
        for readNum in xrange((COVERAGE * len(ref)) / READ_LEN):
        # for readNum in xrange(200):
            start = randint(0, len(ref) - READ_LEN)
            origRead = ref[start:start+READ_LEN]  # + ref[start+READ_LEN:start+READ_LEN+3].lower()
            read = []
            readPos = 0
            refPos = 0
            while readPos < READ_LEN:
                if randint(1, INDEL_ODDS) == 1:
                    length = randint(1, 3)
                    if randint(0, 1):
                        refPos += length
                    else:
                        for i in xrange(length):
                            read.append(choice(BASES))
                            readPos += 1
                elif randint(1, ERR_ODDS) == 1:
                    newBase = choice(BASES)
                    while newBase == ref[start+refPos]:
                        newBase = choice(BASES)
                    read.append(newBase)
                    refPos += 1
                    readPos += 1
                else:
                    read.append(ref[start+refPos])
                    refPos += 1
                    readPos += 1
            read = ''.join(read[:READ_LEN])
            if randint(0, 1):
                read = reverseComplement(read)
                origRead = reverseComplement(origRead)
            fp.write('@' + str(readNum) + '\n')
            fp.write(read + '\n')
            fp.write('+\n')
            fp.write(origRead + '\n')


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
    fixed = 0
    botched = 0
    tried = 0
    ignored = 0
    with open('/playpen/sgreens/fake/corrected2.fa', 'r') as fp:
        for read, _, perfectRead, _, origRead in grouper(fp, 5):
            if origRead == perfectRead:  # read had no error
                if read == perfectRead:
                    ignored += 1
                else:
                    botched += 1
            else:  # read had error
                if read == perfectRead:
                    fixed += 1
                else:
                    tried += 1
    print 'Fixed', fixed, 'reads out of', fixed + tried, 'reads with errors: %.2f%%' % ((fixed*100.)/(fixed+tried))
    print 'Ignored', ignored, 'reads out of', ignored + botched, 'reads without errors: %.2f%%' % ((ignored*100.)/(ignored+botched))
    print 'Altered', fixed+botched, 'reads out of', fixed + botched + ignored + tried,\
        'total reads: %.2f%%' % (((fixed+botched)*100.)/(fixed+botched+ignored+tried))


def main():
    createReads()
    # correction()
    # measurePerformance()


if __name__ == '__main__':
    main()