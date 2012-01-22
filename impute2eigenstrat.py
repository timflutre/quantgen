#!/usr/bin/env python

# Author: Tim Flutre
# Aim: convert a genotype file in the IMPUTE format into the three files
# 'genotype', 'snp' and 'indiv', required by EIGENSOFT
# Note: the first allele in the IMPUTE file is considered as being the reference allele for EIGENSOFT

import sys
import gzip

if len(sys.argv) < 4:
    msg = "ERROR: missing arguments"
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)

imputeFile = sys.argv[1]       # gzipped, header line
eigensoftPrefix = sys.argv[2]  # prefix for the 3 output files
label = sys.argv[3]            # used in the 'sample' output file (eg. Case, Controls, Pop1)

# make file names for EIGENSOFT
genoFile = "%s_eigenstrat_genotypes.txt" % eigensoftPrefix
snpFile = "%s_eigenstrat_snp.txt" % eigensoftPrefix
indivFile = "%s_eigenstrat_indiv.txt" % eigensoftPrefix

# open file handlers
imputeH = gzip.open(imputeFile)
genoH = open(genoFile, "w")
snpH = open(snpFile, "w")
indivH = open(indivFile, "w")

# make the sample file thanks to the header
line = imputeH.readline()
lSamples = line.split()[5:]
nbSamples = len(lSamples)
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
    if countLines % 50000 == 0: print countLines; sys.stdout.flush()
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
            txt += "9"  # missing data
    genoH.write("%s\n" % txt)
print "nb of SNPs: %i" % countLines
    
# close file handlers
imputeH.close()
indivH.close()
genoH.close()
snpH.close()
