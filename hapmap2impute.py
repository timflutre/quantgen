#! /usr/bin/env python

# Author: Tim Flutre
# License: GPL-3
# Aim: convert a genotype file from the HapMap format into the IMPUTE format


import sys
import os
import getopt
import gzip


def help():
    msg = "`%s' converts a genotype file from the HapMap into the IMPUTE format.\n" % os.path.basename(sys.argv[0])
    msg += "\n"
    msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
    msg += "\n"
    msg += "Options:\n"
    msg += "  -h, --help\tdisplay the help and exit\n"
    msg += "  -v, --verbose\tverbosity level (default=1)\n"
    msg += "  -i, --input\tinput pattern with the whole path to the HapMap file (gzipped, with a header line)\n"
    msg += "  -o, --output\toutput file in the IMPUTE format (gzipped, with a header line)\n"
    msg += "  -c, --chr\tchromosome numbers to convert (eg. '1-17-3', all autosomes by default)\n"
    msg += "  -b, --bed\tcoordinate of SNPs in a BED file (eg. output of liftOver, gzipped)\n"
    msg += "  -s, --snp\tfile with a list of SNP identifiers to ignore (eg. if they have conflicting coordinates)\n"
    msg += "\n"
    msg += "Example:\n"
    msg += "  %s -i ~/HMr28/genotypes_CHR_CEU_r28_nr.b36_fwd.txt.gz -o genotypes_allchrs_CEU.impute.gz" % os.path.basename(sys.argv[0])
    print msg; sys.stdout.flush()
    
    
def setParamsFromCmdLine():
    inPattern = ""
    outFile = ""
    lChrs = range(1,23)
    snpCoordFile = ""
    snpIgnoreFile = ""
    verbose = 1
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv:i:o:c:b:s:",
                                   ["help", "verbose=", "input=",
                                    "output=", "chr=", "bed=", "snp="])
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
            inPattern = a
        elif o in ("-o", "--output"):
            outFile = a
        elif o in ("-c", "--chr"):
            lChrs = map(int, a.split("-"))
        elif o in ("-b", "--bed"):
            snpCoordFile = a
        elif o in ("-s", "--snp"):
            snpIgnoreFile = a
    if inPattern == "":
        msg = "ERROR: missing input pattern (-i)"
        sys.stderr.write("%s\n\n" % msg)
        help()
        sys.exit(1)
    if outFile == "":
        msg = "ERROR: missing output file (-o)"
        sys.stderr.write("%s\n\n" % msg)
        help()
        sys.exit(1)
    return inPattern, outFile, lChrs, snpCoordFile, snpIgnoreFile, verbose


def loadFileWithListOfSnpsToIgnore(snpIgnoreFile, verbose):
    lSnpsToIgnore = []
    
    if snpIgnoreFile != "":
        if verbose > 0:
            print "load SNPs to ignore..."
            sys.stdout.flush()
        
        snpIgnoreH = open(snpIgnoreFile)
        while True:
            line = snpIgnoreH.readline()
            if line == "": break
            if line[:-1] not in lSnpsToIgnore:
                lSnpsToIgnore.append(line[:-1])
        snpIgnoreH.close()
        
        if verbose > 0:
            print "nb of SNPs to ignore: %i" % len(lSnpsToIgnore)
            sys.stdout.flush()
            
    return lSnpsToIgnore


def getNewSnpCoordinates(snpCoordFile, lSnpsToIgnore, verbose):
    dSnpId2NewCoord = {}
    
    if snpCoordFile != "":
        if verbose > 0:
            print "load new SNP coordinates..."
            sys.stdout.flush()
            
        snpCoordH = gzip.open(snpCoordFile)
        
        while True:
            line = snpCoordH.readline()
            if line == "": break
            tok = line[:-1].split()
            if tok[3] in lSnpsToIgnore:
                continue
            if dSnpId2NewCoord.has_key(tok[3]):
                msg = "ERROR: SNP '%s' is redundant" % tok[3]
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            dSnpId2NewCoord[tok[3]] = tok[2]
        snpCoordH.close()
        
        if verbose > 0:
            print "nb of SNPs with new coords: %i" % len(dSnpId2NewCoord.keys())
            sys.stdout.flush()
            
    return dSnpId2NewCoord


def convertHapMapToImpute(inPattern, chrNb, lSnpsToIgnore, dSnpId2NewCoord, outH, isFirstFile, verbose):
    hmFile = inPattern.replace("CHR", "chr%i" % chrNb)
    if verbose > 0:
        print "convert file %s..." % hmFile
        sys.stdout.flush()
        
    hmH = gzip.open(hmFile)
    
    # handle the header (write one only for the first file)
    line = hmH.readline()
    if isFirstFile:
        lToks = line[:-1].split()
        txt = "chr id coord a1 a2"
        for i in range(11,len(lToks)):
            txt += " %s" % lToks[i]
        outH.write("%s\n" % txt)
        
    # handle the other lines
    while True:
        line = hmH.readline()
        if line == "": break
        lToks = line[:-1].split()
        snpId = lToks[0]
        if snpId not in lSnpsToIgnore:
            txt = "chr%i" % chrNb
            txt += " %s" % snpId
            if dSnpId2NewCoord.has_key(snpId):
                txt += " %s" % dSnpId2NewCoord[snpId]
            else:
                txt += " %s" % lToks[2]
            a1, a2 = lToks[1].split("/")
            txt += " %s" % a1
            txt += " %s" % a2
            for i in range(11,len(lToks)):
                if lToks[i] == a1+a1:
                    txt += " 1 0 0"
                elif a1 in lToks[i] and a2 in lToks[i]:
                    txt += " 0 1 0"
                elif lToks[i] == a2+a2:
                    txt += " 0 0 1"
                else:
                    txt += " 0 0 0"
            outH.write("%s\n" % txt)
            
    hmH.close()
    
    
def main():
    
    inPattern, outFile, lChrs, snpCoordFile, snpIgnoreFile, verbose = setParamsFromCmdLine()
    
    lSnpsToIgnore = loadFileWithListOfSnpsToIgnore(snpIgnoreFile, verbose)
    
    dSnpId2NewCoord = getNewSnpCoordinates(snpCoordFile, lSnpsToIgnore, verbose)
    
    outH = gzip.open(outFile, "w")
    isFirstFile = True
    for chrNb in lChrs:
        convertHapMapToImpute(inPattern, chrNb, lSnpsToIgnore, dSnpId2NewCoord, outH, isFirstFile, verbose)
        isFirstFile = False
    outH.close()
    
    
if __name__ == "__main__":
    main()
