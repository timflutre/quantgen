#!/usr/bin/env python

import sys
import os
import getopt
import gzip

def usage():
    msg = "usage: %s [options]" % os.path.basename(sys.argv[0])
    msg += "\n-i input SAM file"
    msg += "\n-e max edit distance (default=0)"
    msg += "\n-o output prefix"
    msg += "\n-v verbose (0/default=1/2)"
    msg += "\naim: keep reads matching uniquely and (almost) perfectly"
    msg += "\nreturn reads complying to filter in SAM and BED files"
    sys.stderr.write( "%s\n" % msg )

inFile = ""      # input SAM file
maxEditDist = 0  # max nb of differences (mismatches/indels)
outPrefix = ""   # output prefix
verbose = 1
try:
    opts, args = getopt.getopt( sys.argv[1:], "hi:e:o:v:" )
except getopt.GetoptError, err:
    sys.stderr.write( str(err)+"\n" ); usage(); sys.exit(2)
for o, a in opts:
    if o == "-h":
        usage(); sys.exit()
    elif o == "-i":
        inFile = a
    elif o == "-e":
        maxEditDist = int(a)
    elif o == "-o":
        outPrefix = a
    elif o == "-v":
        verbose = int(a)
if inFile == "" or outPrefix == "":
    usage(); sys.exit(1)

if verbose > 0:
    msg = " ".join( sys.argv )
    sys.stderr.write( "%s\n" % msg )
    msg = "input SAM: %s" % inFile
    msg += "\nmax edit distance: %i" % maxEditDist
    sys.stderr.write( "%s\n" % msg )
inH = open( inFile )
outSamFile = "%s.sam" % outPrefix
outSamH = open( outSamFile, "w" )
outBedFile = "%s.bed.gz" % outPrefix
outBedH = gzip.open( outBedFile, "w" )
#outBedH.write( "chrom chromStart chromEnd name score strand\n" )
if verbose > 0:
    msg = "outputs: '%s' and '%s'" % ( outSamFile, outBedFile )
    sys.stderr.write( "%s\n" % msg )
    sys.stderr.flush()

countReads = 0
countMapped = 0
countUniq = 0
countPerfect = 0
countAsPerfect = 0
countSaved = 0

# for each input line (~each read)
while True:
    line = inH.readline().replace("\n","")
    if line == "": break
    if line[0] == "@": continue
    countReads += 1
    lTok = line.split("\t")

    # if the read was mapped
    if lTok[1] != "4":
        countMapped += 1

        # if the read is unique and perfect
        isUniq = False
        isAsPerfect = False
        type_XT = ""
        nbBestHits_X0 = 1000000000
        nbSubOptHits_X1 = 1000000000
        editDist_NM = 1000000000
        for opt in lTok[11:]:
            if "XT" in opt:
                type_XT = opt.split(":")[2]
            elif "X0" in opt:
                nbBestHits_X0 = int(opt.split(":")[2])
            elif "X1" in opt:
                nbSubOptHits_X1 = int(opt.split(":")[2])
            elif "NM" in opt:
                editDist_NM = int(opt.split(":")[2])
        if type_XT == "U" and nbBestHits_X0 == 1 and nbSubOptHits_X1 == 0:
            isUniq = True
            countUniq += 1
        if isUniq:
            if editDist_NM == 0:
                isAsPerfect = True
                countPerfect += 1
            else:
                if editDist_NM <= maxEditDist:
                    isAsPerfect = True
                    countAsPerfect += 1
                else:
                    if verbose > 1:
                        msg = "high edit distance", lTok[0], lTok[11:]
                        sys.stderr.write( "%s\n" % msg )

        # save it
        if isUniq and isAsPerfect:
            countSaved += 1

            # in SAM format
            outSamH.write( "%s\n" % line )

            # in BED format
            length = len(lTok[9])
            strand = ""
            if lTok[1] == "0":
                strand = "+"
            elif lTok[1] == "16":
                strand = "-"
            bed = "%s\t%i\t%i" % ( lTok[2], int(lTok[3])-1, int(lTok[3])+length )
            bed += "\t%s\t1000\t%s" % ( lTok[0], strand )
            outBedH.write( "%s\n" % bed )

if verbose > 0:
    msg = "nb of reads: %i" % countReads
    msg += "\nnb of mapped reads: %i" % countMapped
    msg += "\nnb of reads mapping uniquely: %i" % countUniq
    msg += "\nnb of reads mapping perfectly: %i" % countPerfect
    msg += "\nnb of reads mapping almost perfectly: %i" % countAsPerfect
    msg += "\nnb of saved reads: %i" % countSaved
    sys.stderr.write( "%s\n" % msg )
inH.close()
outSamH.close()
