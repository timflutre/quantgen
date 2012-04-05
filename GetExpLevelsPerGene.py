#! /usr/bin/env python

## This script sums all read counts from exons per gene and computes their RPKM.
## Author: Timothee Flutre

import sys
import os
import getopt
import scipy


class Gene(object):
    def __init__(self, n="", c="", s=0, e=0):
        self.name = n
        self.chr = c
        self.start = int(s)
        self.end = int(e)
        self.sumExonLengths = self.end - self.start
        self.aReadCounts = None
    def __repr__(self):
        return self.name, self.chr, self.start, self.end, self.sumExonLengths
    def getRpkm(self, aAllReadCounts):
        rpkms = scipy.zeros(len(aAllReadCounts))
        idx = (aAllReadCounts != 0)
        rpkms[idx] = 10**9 * self.aReadCounts[idx] / (aAllReadCounts[idx] \
                                                          * self.sumExonLengths)
        return rpkms
    def read(self, inH):
        line = inH.readline()
        if line == "":
            return
        lTokens = line.rstrip().split()
        self.name = lTokens[3]
        self.chr = lTokens[0]
        self.start = int(lTokens[1])
        self.end = int(lTokens[2])
        self.sumExonLengths = self.end - self.start
        self.aReadCounts = scipy.array(map(int, lTokens[4:]))
        while True:
            prevPos = inH.tell()
            line = inH.readline()
            if line == "":
                break
            lTokens = line.rstrip().split()
            if lTokens[3] != self.name:
                inH.seek(prevPos)
                break
            if int(lTokens[1]) < self.start:
                self.start = int(lTokens[1])
            if int(lTokens[2]) > self.end:
                self.end = int(lTokens[2])
            self.sumExonLengths += int(lTokens[2]) - int(lTokens[1])
            self.aReadCounts += scipy.array(map(int, lTokens[4:]))
    def write(self, outH, expLevels, aAllReadReadCounts):
        txt = "%s" % self.chr
        txt += " %i" % self.start
        txt += " %i" % self.end
        txt += " %s" % self.name
        if expLevels == "rpkm":
            rpkms = self.getRpkm(aAllReadReadCounts)
            for i in xrange(0,len(rpkms)):
                txt += " %.2f" % rpkms[i]
        else:
            for i in xrange(0,len(self.aReadCounts)):
                txt += " %i" % self.aReadCounts[i]
        outH.write("%s\n" % txt)
        
        
def help():
    msg = "`%s' sums all read counts from exons per gene and computes their RPKM.\n" % os.path.basename(sys.argv[0])
    msg += "RPKM = 10^9 * C / ( N * L )\n"
    msg += "C=read count for transcript, N=total read count, L=sum of exon lengths in bp\n"
    msg += "\n"
    msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
    msg += "\n"
    msg += "Options:\n"
    msg += " -h, --help\tdisplay the help and exit\n"
    msg += " -v, --verbose\tverbosity level (default=1)\n"
    msg += " -i, --input\tinput file in BED-like format\n"
    msg += "\t\tcolumns: chr exonStart exonEnd geneName readCountSample1, readCountSample2...\n"
    msg += "\t\tfile should be sorted, use sort -k1,1V -k4,4V -k2,2n -k3,3n\n"
    msg += "\t\tshould have a header line with the identifiers of the samples\n"
    msg += " -o, --output\toutput file in a similar format\n"
    msg += "\t\tcolumns: chr most5'exon most3'exon geneName RpkmSample1 RpkmSample2...\n"
    msg += " -e, --exp\texpression levels as 'rpkm' (default) or 'readcount'\n"
    msg += "\n"
    msg += "Example:\n"
    msg += " %s -i read_counts_per_exon.txt -o rpkm_per_gene.txt" % os.path.basename(sys.argv[0])
    print msg; sys.stdout.flush()
    
    
def setParamsFromCmdLine():
    inFile = ""
    outFile = ""
    expLevels = "rpkm"
    verbose = 1
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv:i:o:e:",
                                   ["help", "verbose=", "input=",
                                    "output=", "exp="])
    except getopt.GetoptError, err:
        sys.stderr.write("%s\n" % str(err))
        help()
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            help()
            sys.exit(0)
        elif o in ("-v", "--verbose"):
            verbose = int(a)
        elif o in ("-i", "--input"):
            inFile = a
        elif o in ("-o", "--output"):
            outFile = a
        elif o in ("-e", "--exp"):
            expLevels = a
    if inFile == "":
        msg = "ERROR: missing input file (-i)"
        sys.stderr.write("%s\n\n" % msg)
        help()
        sys.exit(1)
    if outFile == "":
        msg = "ERROR: missing output file (-o)"
        sys.stderr.write("%s\n\n" % msg)
        help()
        sys.exit(1)
    if expLevels not in ("rpkm", "readcount"):
        msg = "ERROR: unknown expression level (-e)"
        sys.stderr.write("%s\n\n" % msg)
        help()
        sys.exit(1)
    return inFile, outFile, expLevels, verbose
    
    
def readAllGenes(inFile, verbose):
    if verbose > 0:
        print "load file '%s'..." % inFile
        sys.stdout.flush()
    lSampleNames = []
    lGenes = []
    aAllReadCounts = None
    inH = open(inFile)
    line = inH.readline()
    lTokens = line.rstrip().split()
    lSampleNames = lTokens[4:]
    aAllReadCounts = scipy.zeros(len(lSampleNames))
    lGeneNames = []
    while True:
        iGene = Gene()
        iGene.read(inH)
        if iGene.name == "":
            break
        if iGene.name in lGeneNames:
            msg = "ERROR: the input file is not sorted, see the help"
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
        aAllReadCounts += iGene.aReadCounts
        lGenes.append(iGene)
        lGeneNames.append(iGene.name)
    inH.close()
    if verbose > 0:
        print "nb of genes: %i" % len(lGenes)
        print "nb of samples: %i" % len(aAllReadCounts)
        print "nb of reads: %i" % aAllReadCounts.sum()
        sys.stdout.flush()
    return lSampleNames, lGenes, aAllReadCounts


def writeAllGenes(lSampleNames, lGenes, aAllReadCounts, expLevels, outFile, verbose):
    outH = open(outFile, "w")
    txt = "chr start end gene"
    for sample in lSampleNames:
        txt += " %s" % sample
    outH.write("%s\n" % txt)
    for iGene in lGenes:
        iGene.write(outH, expLevels, aAllReadCounts)
    outH.close()
    
    
def main():
    inFile, outFile, expLevels, verbose = setParamsFromCmdLine()
    lSampleNames, dGenes, aAllReadCounts = readAllGenes(inFile, verbose)
    writeAllGenes(lSampleNames, dGenes, aAllReadCounts, expLevels, outFile, verbose)
    
    
if __name__ == "__main__":
    main()
