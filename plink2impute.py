#!/usr/bin/env python

# Author: Timothee Flutre
# License: GPL-3
# Aim: converts PLINK genotype data into the IMPUTE format
# help2man -o MyClass.man ./MyClass.py
# groff -mandoc MyClass.man > MyClass.ps


import sys
import os
import getopt
import time
import datetime
import math
import gzip


class Plink2Impute(object):
    
    def __init__(self):
        self.verbose = 1
        self.tpedFile = ""
        self.tfamFile = ""
        self.rmvMarkerFile = ""
        self.inSampleFile = ""
        self.newCoordFile = ""
        self.keepOnlyAutosomes = True
        self.outPrefix = "genotypes"
        self.outSuffix = ".study.gens"
        
        
    def help(self):
        msg = "`%s' converts PLINK genotype data into the IMPUTE format.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += " -h, --help\tdisplay the help and exit\n"
        msg += " -V, --version\toutput version information and exit\n"
        msg += " -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += " -t\t\tfile with the genotypes (in TPED format)\n"
        msg += " -m\t\tfile with the samples (in TFAM format)\n"
        msg += " -s\t\tfile with samples to keep\n"
        msg += "\t\tsamples will be saved in their order in this file\n"
        msg += "\t\tif two columns (pop<tab>sample), outputs will be 'pop...'\n"
        msg += " -c\t\tfile with new marker coordinates (in BED format)\n"
        msg += "\t\teg. output from liftOver\n"
        msg += " -r\t\tfile with a list of markers to remove\n"
        msg += "\t\tone marker identifier per line\n"
        msg += "\t\tlines containing '#' are skipped\n"
        msg += " -k\t\tkeep all chromosomes (and not only autosomes)\n"
        msg += " -o\t\toutput suffix after 'chr?' (default='.study.gens')\n"
        msg += "\t\toutput prefix is always 'genotypes'\n"
        msg += "\t\toutput files have a header line and are not gzipped\n"
        msg += "\n"
        msg += "Remarks:\n"
        msg += " If -c but not -r, markers with no new coordinates will be skipped.\n"
        msg += " Handle gzipped files as long as they finish in 'gz'.\n"
        print msg; sys.stdout.flush()
        
        
    def version(self):
        msg = "%s 0.1\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "This is free software; see the source for copying conditions.  There is NO\n"
        msg += "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
        print msg; sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], "hVv:t:m:s:c:r:ko:",
                                       ["help", "version", "verbose="])
        except getopt.GetoptError, err:
            sys.stderr.write("%s\n" % str(err))
            self.help()
            sys.exit(2)
        for o, a in opts:
            if o == "-h" or o == "--help":
                self.help()
                sys.exit(0)
            elif o == "-V" or o == "--version":
                self.version()
                sys.exit(0)
            elif o == "-v" or o == "--verbose":
                self.verbose = int(a)
            elif o == "-t":
                self.tpedFile = a
            elif o == "-m":
                self.tfamFile = a
            elif o == "-r":
                self.rmvMarkerFile = a
            elif o == "-s":
                self.inSampleFile = a
            elif o == "-c":
                self.newCoordFile = a
            elif o == "-k":
                self.keepOnlyAutosomes = False
            elif o == "-o":
                self.outSuffix = a
            else:
                assert False, "unhandled option"
                
                
    def checkAttributes(self):
        if self.tpedFile == "":
            msg = "ERROR: missing required argument -t"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.tfamFile == "":
            msg = "ERROR: missing required argument -m"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not self.keepOnlyAutosomes:
            msg = "ERROR: keeping sex chromosomes not yet available"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
            
            
    def getMarkersToRemove(self):
        lMarkersToRmv = []
        if self.rmvMarkerFile != "":
            if self.verbose > 0:
                print "get markers to remove ..."
                sys.stdout.flush()
            rmvMarkerH = open(self.rmvMarkerFile)
            while True:
                line = rmvMarkerH.readline().rstrip()
                if line == "":
                    break
                if "#" in line:
                    continue
                if line not in lMarkersToRmv:
                    lMarkersToRmv.append(line)
            rmvMarkerH.close()
            if self.verbose > 0:
                print "nb of markers to remove: %i" % len(lMarkersToRmv)
                sys.stdout.flush()
        return lMarkersToRmv
    
    
    def getNewMarkerCoordinates(self):
        dMarker2NewCoord = {}
        
        if self.newCoordFile != "":
            if self.verbose > 0:
                print "get new marker coordinates ..."
                sys.stdout.flush()
            if self.newCoordFile[-2:] == "gz":
                newCoordH = gzip.open(self.newCoordFile)
            else:
                newCoordH = open(self.newCoordFile)
                
            while True:
                line = newCoordH.readline()
                if line == "":
                    break
                if "\t" in line:
                    tokens = line.split("\t")
                else:
                    tokens = line.split(" ")
                if dMarker2NewCoord.has_key(tokens[3]):
                    msg = "ERROR: marker '%s' is redundant" % tokens[3]
                    sys.stderr.write("%s\n" % msg)
                    sys.exit(1)
                dMarker2NewCoord[tokens[3]] = tokens[2]
                
            newCoordH.close()
            if self.verbose > 0:
                print "nb of markers with new coordinates: %i" % \
                    len(dMarker2NewCoord)
                sys.stdout.flush()
                
        return dMarker2NewCoord
    
    
    def getSamplesToKeep(self):
        lSamplesToKeep = []
        dPop2Samples = {}
        
        if self.inSampleFile != "":
            if self.verbose > 0:
                print "get samples to keep ..."
                sys.stdout.flush()
            inSampleH = None
            if self.inSampleFile[-2:] == "gz":
                inSampleH = gzip.open(self.inSampleFile)
            else:
                inSampleH = open(self.inSampleFile)
                
            while True:
                line = inSampleH.readline().rstrip()
                if line == "":
                    break
                if line in lSamplesToKeep:
                    msg = "ERROR: sample '%s' is redundant" % line
                    sys.stderr.write("%s\n" % msg)
                    sys.exit(1)
                tokens = line.split()
                lSamplesToKeep.append(tokens[0])
                if len(tokens) > 1:
                    if not dPop2Samples.has_key(tokens[1]):
                        dPop2Samples[tokens[1]] = []
                    dPop2Samples[tokens[1]].append(tokens[0])
                    
            inSampleH.close()
            if self.verbose > 0:
                print "nb of samples to keep: %i" % len(lSamplesToKeep)
                if dPop2Samples != {}:
                    print "nb of populations: %i" % len(dPop2Samples)
                sys.stdout.flush()
                
        return lSamplesToKeep, dPop2Samples
    
    
    def getSampleInitialOrderFromTfamFile(self, lSamplesToKeep):
        dSample2LineIdx = {}
        if self.tfamFile[-2:] == "gz":
            tfamH = gzip.open(self.tfamFile)
        else:
            tfamH = open(self.tfamFile)
            
        lineIdx = 0
        while True:
            line = tfamH.readline().rstrip()
            if line == "":
                break
            tokens = line.split(" ")
            if tokens[1] in lSamplesToKeep:
                dSample2LineIdx[tokens[1]] = lineIdx
            lineIdx += 1
            
        if len(dSample2LineIdx) != len(lSamplesToKeep):
            msg = "ERROR: some samples to keep are absent from the TFAM file"
            sys.stderr.write("%s\n" % msg)
            for sample in lSamplesToKeep:
                if not dSample2LineIdx.has_key(sample):
                    print sample
            sys.exit(1)
            
        tfamH.close()
        return dSample2LineIdx
    
    
    def convertTpedFile(self, lMarkersToRmv, dMarker2NewCoord,
                        lSamplesToKeep, dPop2Samples, dSample2LineIdx):
        if self.verbose > 0:
            print "convert TPED file ..."
            sys.stdout.flush()
            
        outH = None
        dPop2OutH = {}
        txt = "chr id coord a1 a2"
        if dPop2Samples == {}:
            outFile = "%s%s" % (self.outPrefix, self.outSuffix)
            outH = open(outFile, "w")
        else:
            for pop in dPop2Samples:
                outFile = "%s_%s%s" % (self.outPrefix, pop, self.outSuffix)
                outH = open(outFile, "w")
                dPop2OutH[pop] = outH
                dPop2OutH[pop].write("%s" % txt)
                txtSamples = ""
                for sample in dPop2Samples[pop]:
                    txtSamples += " %s" % sample
                dPop2OutH[pop].write("%s\n" % txtSamples)
                
        if self.tpedFile[-2:] == "gz":
            tpedH = gzip.open(self.tpedFile)
        else:
            tpedH = open(self.tpedFile)
            
        nbMarkers = 0
        nbMarkersSaved = 0
        
        while True:
            line = tpedH.readline().rstrip()
            if line == "":
                break
            
            nbMarkers += 1
#            if nbMarkers == 3: break
            if self.verbose > 0 and nbMarkers % 100000 == 0:
                print nbMarkers; sys.stdout.flush()
                
            tokens = line.split()
            
            chrNb = int(tokens[0])
            if self.keepOnlyAutosomes and chrNb > 22:
                if self.verbose > 0:
                    print "skip marker '%s' because not on autosome" % marker
                    sys.stdout.flush()
                continue
            marker = tokens[1]
            if marker in lMarkersToRmv:
                if self.verbose > 0:
                    print "skip marker '%s' because flagged as removed" % marker
                    sys.stdout.flush()
                continue
            if lMarkersToRmv == [] and not dMarker2NewCoord.has_key(marker):
                if self.verbose > 0:
                    print "skip marker '%s' because no new coordinate" % marker
                    sys.stdout.flush()
                continue
            nbMarkersSaved += 1
            
            if dMarker2NewCoord.has_key(marker):
                coord = dMarker2NewCoord[marker]
            else:
                coord = tokens[3]
            txt = "chr%i %s %s" % (chrNb, marker, coord)
            
            sAlleles = set(tokens[4:])
            sAlleles.discard("0")
            sAlleles.discard("N")
            if len(sAlleles) == 0:
                msg = "ERROR: no allele for marker '%s'" % marker
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            if len(sAlleles) == 1:
                msg = "ERROR: no polymorphism for marker '%s'" % marker
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            elif len(sAlleles) == 2:
                A = sAlleles.pop()
                B = sAlleles.pop()
            else:
                msg = "ERROR: too many alleles for '%s'" % marker
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            txt += " %s %s" % (A, B)
            
            if dPop2Samples == {}:
                outH.write("%s" % txt)
                txtGeno = ""
                for sample in lSamplesToKeep:
                    colIdx = 4 + 2 * dSample2LineIdx[sample]
                    if tokens[colIdx] + tokens[colIdx+1] == A + A:
                        txtGeno += " 1 0 0"
                    elif tokens[colIdx] + tokens[colIdx+1] == A + B \
                            or tokens[colIdx] + tokens[colIdx+1] == B + A:
                        txtGeno += " 0 1 0"
                    else:
                        txtGeno += " 0 0 1"
                outH.write("%s\n" % txtGeno)
            else:
                for pop, outH in dPop2OutH.iteritems():
                    outH.write("%s" % txt)
                    txtGeno = ""
                    for sample in dPop2Samples[pop]:
                        colIdx = 4 + 2 * dSample2LineIdx[sample]
                        if tokens[colIdx] + tokens[colIdx+1] == A + A:
                            txtGeno += " 1 0 0"
                        elif tokens[colIdx] + tokens[colIdx+1] == A + B \
                                or tokens[colIdx] + tokens[colIdx+1] == B + A:
                            txtGeno += " 0 1 0"
                        else:
                            txtGeno += " 0 0 1"
                    outH.write("%s\n" % txtGeno)
                    
        tpedH.close()
        for pop in dPop2OutH:
            dPop2OutH[pop].close()
            
        if self.verbose > 0:
            print "nb of input markers: %i" % nbMarkers
            print "nb of output markers: %i" % nbMarkersSaved
            sys.stdout.flush()
            
            
    def run(self):
        self.checkAttributes()
        
        if self.verbose > 0:
            msg = "START %s" % time.strftime("%Y-%m-%d %H:%M:%S")
            startTime = time.time()
            print msg; sys.stdout.flush()
            
        lMarkersToRmv = self.getMarkersToRemove()
        
        dMarker2NewCoord = self.getNewMarkerCoordinates()
        
        lSamplesToKeep, dPop2Samples = self.getSamplesToKeep()
        
        dSample2LineIdx = self.getSampleInitialOrderFromTfamFile(lSamplesToKeep)
        
        self.convertTpedFile(lMarkersToRmv, dMarker2NewCoord,
                             lSamplesToKeep, dPop2Samples,
                             dSample2LineIdx)
        
        if self.verbose > 0:
            msg = "END %s" % time.strftime("%Y-%m-%d %H:%M:%S")
            endTime = time.time()
            runLength = datetime.timedelta(seconds=
                                           math.floor(endTime - startTime))
            msg += " (%s)" % str(runLength)
            print msg; sys.stdout.flush()
            
            
if __name__ == "__main__":
    i = Plink2Impute()
    i.setAttributesFromCmdLine()
    i.run()
