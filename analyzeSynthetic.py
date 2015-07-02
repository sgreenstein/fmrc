from MUSCython import MultiStringBWTCython as msbwt
from itertools import izip
from collections import Counter
from string import maketrans

TRAN_TAB = maketrans('ACGT', 'TGCA')


def grouper(iterable, n):
    args = [iter(iterable)] * n
    return izip(*args)


def revComp(read):
    return read[::-1].translate(TRAN_TAB)


def counts(bwt, read, k):
    counts = []
    for start, end in zip(xrange(len(read) - k), xrange(k, len(read))):
        lo, hi = bwt.findIndicesOfStr(read[start:end])
        counts.append(hi - lo)
    return counts


def main():
    bwt = msbwt.loadBWT('/playpen/sgreens/ecoli/msbwtfake/bwt_50')
    crcSubs = Counter()
    missedSubs = Counter()
    falsePos = 0
    trueNeg = 0
    crcDels = 0
    missedDels = 0
    crcIns = 0
    missedIns = 0
    with open('/playpen/sgreens/ecoli/msbwtfake/corrected.fastq') as fp:
    # with open('/playpen/sgreens/ecoli/blessfake/corrected.corrected.fastq') as fp:
        for readName, read, errString, origRead in grouper(fp, 4):
            read = read[:-1]
            origRead = origRead[:-1]
            errString = errString[:-1]
            numSubs = errString.count('S')
            numDels = errString.count('D')
            numIns = errString.count('I')
            if read not in origRead:
                for base, origBase, err in zip(read, origRead[3:-3], errString):
                    if err == 'D':
                        firstIndex = errString.find('D')
                        lastIndex = firstIndex + errString[firstIndex:].replace('I', '-').replace('S', '-').find('-')
                        if read[firstIndex:lastIndex] == origRead[3:-3][firstIndex:lastIndex]:
                            crcDels += 1
                        else:
                            missedDels += 1
                        break
                    if err == 'I':
                        firstIndex = errString.find('I')
                        lastIndex = firstIndex + errString[firstIndex:].replace('D', '-').replace('S', '-').find('-')
                        if read[firstIndex:lastIndex] == origRead[3:-3][firstIndex:lastIndex]:
                            crcIns += 1
                        else:
                            missedIns += 1
                        break
                    elif err == 'S':
                        if base == origBase:
                            crcSubs[numSubs] += 1
                        else:
                            missedSubs[numSubs] += 1
                    elif err == '-':
                        if base != origBase:
                            falsePos += 1
                        else:
                            trueNeg += 1
            else:
                if numDels > 0:
                    crcDels += 1
                if numIns > 0:
                    crcIns += 1
                if numSubs > 0:
                    crcSubs[numSubs] += numSubs
                trueNeg += errString.count('-')
    totalCorrected = 0
    for numSubs, corrected in crcSubs.iteritems():
        totalCorrected += corrected
        print 'In reads with', numSubs, 'subs, corrected',
        print (100. * corrected) / (corrected + missedSubs[numSubs]), '%'
    totalSubs = totalCorrected
    for numSubs, missed in missedSubs.iteritems():
        totalSubs += missed
    print 'Overall corrected', (100. * totalCorrected) / totalSubs, '% of substitutions'
    print 'Specificity (rate of classifying non-errors correctly)', (100. * trueNeg) / (trueNeg + falsePos), '%'
    print 'Corrected', (100. * crcDels) / (crcDels + missedDels), '% of deletions'
    print 'Corrected', (100. * crcIns) / (crcIns + missedIns), '% of insertions'
    print 'Overall corrected',\
        (100. * (totalCorrected + crcDels + crcIns)) / (totalSubs + missedDels + missedIns + crcDels + crcIns),\
        '% of errors'


if __name__ == '__main__':
    main()
