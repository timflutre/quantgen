#!/usr/bin/env python

# Author: Tim Flutre
# License: GPL-3
# Aim: convert a genotype file in the IMPUTE format into the three files
# 'genotype', 'snp' and 'indiv', required by EIGENSOFT
# Note: the first allele in the IMPUTE file is considered as being the reference allele for EIGENSOFT

import sys
import gzip
import getopt
import os


def help():
    msg = "`%s' converts a genotype file in the IMPUTE format\n" % os.path.basename(sys.argv[0])
    msg += "into the three files required by EIGENSOFT.\n"
    msg += "\n"
    msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
    msg += "\n"
    msg += "Options:\n"
    msg += "  -h, --help\tdisplay the help and exit\n"
    msg += "  -v, --verbose\tverbosity level (default=1)\n"
    msg += "  -i, --input\tinput file in the IMPUTE format (gzipped, with a header line)\n"
    msg += "\t\tsee below the description of the header line\n"
    msg += "  -o, --output\tprefix of the output files\n"
    msg += "\t\twill be followed by '_eigenstrat_genotypes/snp/indiv.txt'\n"
    msg += "  -l, --label\tlabel used in the third column of the 'indiv' file\n"
    msg += "\n"
    msg += "Remarks:\n"
    msg += "  The first allele in the IMPUTE file is considered as being the reference\n"
    msg += "allele for EIGENSOF.\n"
    msg += "  The header of the IMPUTE file should contain, in this order, the chromosome\n"
    msg += "identifier, the SNP identifier, its genomic coordinate, the allele coded A,\n"
    msg += "the allele coded B, then followed by the samples identifier. Note that the header\n"
    msg += "of columns in the header is smaller than the number of columns in the rest of the file."
    print msg; sys.stdout.flush()
    
    
def setParamsFromCmdLine():
    verbose = 1
    imputeFile = ""
    eigensoftPrefix = ""
    label = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv:i:o:l:",
                                   ["help", "verbose=", "input=",
                                    "output=", "label="])
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
            imputeFile = a
        elif o in ("-o", "--output"):
            eigensoftPrefix = a
        elif o in ("-l", "--label"):
            label = a
    return verbose, imputeFile, eigensoftPrefix, label


# get command-line arguments
verbose, imputeFile, eigensoftPrefix, label = setParamsFromCmdLine()
if imputeFile == "":
    msg = "ERROR: missing input file (-i)"
    sys.stderr.write("%s\n\n" % msg)
    help()
    sys.exit(1)
if eigensoftPrefix == "":
    msg = "ERROR: missing output prefix (-o)"
    sys.stderr.write("%s\n\n" % msg)
    help()
    sys.exit(1)
if label == "":
    msg = "ERROR: missing label (-l)"
    sys.stderr.write("%s\n\n" % msg)
    help()
    sys.exit(1)
    
# make file names for EIGENSOFT
genoFile = "%s_eigenstrat_genotypes.txt" % eigensoftPrefix
snpFile = "%s_eigenstrat_snp.txt" % eigensoftPrefix
indivFile = "%s_eigenstrat_indiv.txt" % eigensoftPrefix

# open file handlers
imputeH = gzip.open(imputeFile)
genoH = open(genoFile, "w")
snpH = open(snpFile, "w")
indivH = open(indivFile, "w")

# make the indiv file thanks to the header
line = imputeH.readline()
lSamples = line.split()[5:]
nbSamples = len(lSamples)
if verbose > 0:
    print "nb of samples: %i" % nbSamples
for sample in lSamples:
    if len(sample) > 39:
        msg = "ERROR: sample name should be less than 39 characters"
        sys.stderr.write("%s\n" % msg)
        sys.exit(1)
    txt = "%s\tU\t%s" % (sample, label)
    indivH.write("%s\n" % txt)
    
# convert the rest of the data
countLines = 1
while True:
    line = imputeH.readline()
    if line == "": break
    countLines += 1
    if verbose > 0:
        if countLines % 50000 == 0:
            print countLines; sys.stdout.flush()
    tok = line.split(" ")
    
    # write the SNP file
    snpName = tok[1]
    chrId = tok[0].replace("chr","")
    if chrId == "X":
        chrId = "23"
    elif chrId == "Y":
        chrId = "24"
    elif chrId == "mtDNA":
        chrId = "90"
    elif chrId == "XY" or chrId == "YX":
        chrId = "91"
    phyPos = tok[2]
    refAllele = tok[3]
    varAllele = tok[4]
    txt = "%s\t%s\t0.0\t%s\t%s\t%s" % (snpName, chrId, phyPos, refAllele, varAllele)
    snpH.write("%s\n" % txt)
    
    # write the geno file
    txt = ""
    for i in range(0,nbSamples):
        if tok[5+3*i] == "1" and tok[5+3*i+1] == "0" and tok[5+3*i+2] == "0":
            txt += "2"
        elif tok[5+3*i] == "0" and tok[5+3*i+1] == "1" and tok[5+3*i+2] == "0":
            txt += "1"
        elif tok[5+3*i] == "0" and tok[5+3*i+1] == "0" and tok[5+3*i+2] == "1":
            txt += "0"
        else:
            txt += "9"  # missing data -> EIGENSOFT doesn't handle imputed data for the moment
    genoH.write("%s\n" % txt)
if verbose > 0:
    print "nb of SNPs: %i" % countLines
    
# close file handlers
imputeH.close()
indivH.close()
genoH.close()
snpH.close()
