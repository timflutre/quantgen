#!/usr/bin/env python

# Author: Tim Flutre
# License: GPL-3
# Aim: launch IMPUTE2 on a whole-chromosome (split into chunks, and concatenate the results)
# help2man -o LaunchImpute2.man ./LaunchImpute2.py


import sys
import os
import getopt
import time
import datetime
import math
import glob


class LaunchImpute2(object):
    
    def __init__(self):
        self.refDir = ""
        self.lRefFileStemPatterns = []
        self.gmapDir = ""
        self.studyDir = ""
        self.fixStrand = False
        self.Ne = 20000
        self.chunkLength = 5000000  # 5Mb
        self.force = False
        self.gzipOut = False
        self.clean = False
        self.verbose = 1
        self.debug = False
        self.lChrNbs = range(1,23)
        self.templateGmapFile = ""
        
        
    def help(self):
        msg = "`%s' launches IMPUTE v2 on a whole-chromosome\n" % os.path.basename(sys.argv[0])
        msg += "by splitting the data into chunks, and concatenating the results.\n"
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += " -h, --help\tdisplay the help and exit\n"
        msg += " -V, --version\toutput version information and exit\n"
        msg += " -v, --verbose\tverbosity level (default=1)\n"
        msg += " -r\t\tpath to the directory with the reference data\n"
        msg += "\t\te.g. '~/data/hapmap3_r2_plus_1000g_jun2010_b36_ceu/'\n"
        msg += "\t\tor '~/data/1000Genomes.Dec2010.haplotypes_b37/'\n"
        msg += "\t\tor '~/data/ALL_1000G_phase1interim_jun2011_impute/'\n"
        msg += " --rfsp\t\treference file stem pattern (in the dir given by '-r')\n"
        msg += "\t\te.g. 'AFR.CHR.impute.hap'\n"
        msg += "\t\tor 'pilot1.jun2010.b36.CEU.CHR.snpfilt.haps'\n"
        msg += "\t\tor 'ALL_1000G_phase1interim_jun2011_CHR_impute.hap.gz'\n"
        msg += "\t\tor '<rfsp1>,<rfsp2>'\n"
        msg += "\t\tthe '.legend' files should also be there\n"
        msg += " -g\t\tpath to the directory with the genetic maps\n"
        msg += "\t\te.g. 'genetic_maps_b37/'\n"
        msg += "\t\tdefault is directory given by '-r'\n"
        msg += " -s\t\tpath to the directory with the study data\n"
        msg += "\t\tdefault is current directory\n"
        msg += "\t\tshould contain files named 'chr1.study.gens', ...\n"
        msg += " --fix-strand\tstrand argument for IMPUTE2\n"
        msg += "\t\t'-fix_strand_g', otherwise assume files exist\n"
        msg += "\t\te.g. '-strand_g chr1.study.strand' and so on\n"
        msg += " --Ne\t\teffective population size (default=%i)\n" % self.Ne
        msg += " -l\t\tchunk length (default=%i)\n" % self.chunkLength
        msg += " -c\t\tchromosome (e.g. '12', all autosomes by default)\n"
        msg += " --force\tforce to recompute each chunk\n"
        msg += " --oz\t\twrite output as gzipped files\n"
        msg += " --clean\tclean (remove intermediary chunk files)\n"
        msg += " --debug\tfor debugging purposes\n"
        msg += "\t\trun only the two first chunks of each chr\n"
        msg += "Examples:\n"
        msg += "for chrNb in {1..22}; do echo \"LaunchImpute2.py -r 1000Genomes.Dec2010.haplotypes_b37/ --rfsp EUR.CHR.impute.hap -g genetic_maps_b37/ -s CEU_b37_impute2/ --force -c ${chrNb} -v 1\" | qsub -cwd -j y -V -l h_vmem=10g -N \"job_chr\"$chrNb; done\n"
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
            opts, args = getopt.getopt( sys.argv[1:], "hVv:r:g:s:l:c:",
                                        ["help", "version", "verbose=", 
                                         "rfsp=", "fix-strand", "Ne=", 
                                         "force", "oz", "clean", "debug"] )
        except getopt.GetoptError, err:
            sys.stderr.write("%s\n" % str(err))
            self.help()
            sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                self.help()
                sys.exit(0)
            elif o in ("-V", "--version"):
                self.version()
                sys.exit(0)
            elif o in ("-v", "--verbose"):
                self.verbose = int(a)
            elif o == "-r":
                self.refDir = a
            elif o == "--rfsp":
                self.lRefFileStemPatterns = a.split(",")
            elif o == "-g":
                self.gmapDir = a
            elif o == "-s":
                self.studyDir = a
            elif o == "--fix-strand":
                self.fixStrand = True
            elif o == "--Ne":
                self.Ne = int(a)
            elif o == "-l":
                self.chunkLength = int(a)
            elif o == "-c":
                self.lChrNbs = [int(a)]
            elif o == "--force":
                self.force = True
            elif o == "--oz":
                self.gzipOut = True
            elif o == "--clean":
                self.clean = True
            elif o == "--debug":
                self.debug = True
            else:
                assert False, "unhandled option"
                
                
    def checkAttributes(self):
        if self.refDir == "":
            msg = "ERROR: missing reference directory"
            sys.stderr.write("%s\n" % msg)
            self.help()
            sys.exit(1)
        if self.studyDir == "":
            self.studyDir = os.getcwd()
        if self.gmapDir == "":
            self.gmapDir = self.refDir
        if self.verbose > 0:
            print "reference directory: %s" % self.refDir
            print "study directory: %s" % self.studyDir
            print "genetic maps directory: %s" % self.gmapDir
            if self.fixStrand:
                print "strand argument: -fix_strand_g"
            print "Ne: %i" % self.Ne
            print "chunk length: %i" % self.chunkLength
            if self.lChrNbs == []:
                print "chromosomes: all"
            else:
                print "chromosome: %i" % self.lChrNbs[0]
            print "force: %s" % self.force
            print "gzipped output: %s" % self.gzipOut
            print "clean: %s" % self.clean
            print "verbose: %i" % self.verbose
            print "debug: %s" % self.debug
            print "host: %s" % os.uname()[1]
            sys.stdout.flush()
            
            
    def getTemplateFiles(self):
        lGmapFiles = glob.glob("%s/genetic_map*.txt" % self.gmapDir)
        if len(lGmapFiles) == 0:
            msg = "ERROR: can't find genetic map files"
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
        if "b36" in lGmapFiles[0]:
            self.templateGmapFile = "genetic_map_CHR_combined_b36.txt"
        elif "b37" in lGmapFiles[0]:
            self.templateGmapFile = "genetic_map_CHR_combined_b37.txt"
        else:
            msg = "ERROR: check the names of the genetic map files"
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
            
            
    def getListCoordRangesForEachChunk(self, chrName, chunkLength):
        lCoordRanges = []
        
        studyFile = "%s/%s.study.gens" % (self.studyDir, chrName)
        studyH = open(studyFile)
        line = studyH.readline()
        startCoord = int(line.split(" ")[2])
        endCoord = startCoord + chunkLength
        prevCoord = startCoord
        while True:
            line = studyH.readline()
            if line == "": break
            coord = int(line.split(" ")[2])
            if coord < prevCoord:
                msg = "ERROR: input file '%s' is not sorted" % studyFile
                msg += "\nuse 'for i in {1..22}; do sort -t \" \" -k3,3n chr${i}.study.gens -o chr${i}.study.gens; done'"
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            elif coord == prevCoord:
                msg = "warning: '%s' with redundant coordinate '%i', skip it" % (line.split(" ")[0], coord)
                sys.stderr.write("%s\n" % msg)
                sys.stderr.flush()
            elif coord > endCoord:
                lCoordRanges.append( [ startCoord, endCoord ] )
                startCoord = coord
                endCoord = startCoord + chunkLength
            prevCoord = coord
        if lCoordRanges[-1] != [ startCoord, endCoord ]:
            lCoordRanges.append( [ startCoord, endCoord ] )
        studyH.close()

        if self.verbose > 1 and self.debug:
            print "chunks:", lCoordRanges
            
        return lCoordRanges
    
    
    def launchChunkImputation(self, chrName, chunkId, startCoord, endCoord):
        """
        Launch IMPUTE2 on a given chunk.
        """
        outF = "%s.chunk%i" % (chrName, chunkId)
        outF += ("_%i_%i.impute2" % (startCoord, endCoord)).replace("+","")
        if self.verbose > 0:
            msg = "chunk %i" % (chunkId)
            msg += " (%i->%i)..." % (startCoord, endCoord)
            sys.stdout.write(msg)
            sys.stdout.flush()
        if not self.force and os.path.exists( outF ):
            msg = " already imputed\n"
            sys.stdout.write(msg)
            return outF
        beginTime = time.time()
        tmpF = "impute2.current_%s" % chrName
        cmd = "impute2"
        cmd += " -m %s/" % self.gmapDir
        cmd += self.templateGmapFile.replace("CHR",chrName)

        cmd += " -h"
        for rfsp in self.lRefFileStemPatterns:
            cmd += " %s/" % self.refDir
            cmd += rfsp.replace("CHR",chrName) 
        cmd += " -l"
        for rfsp in self.lRefFileStemPatterns:
            cmd += " %s/" % self.refDir
            prefix = ".".join(rfsp.replace("CHR",chrName).split(".hap")[:-1])
            cmd += "%s.legend" % prefix
            
        cmd += " -g %s/%s.study.gens" % (self.studyDir, chrName)
        if self.fixStrand:
            cmd += " -fix_strand_g"
        else:
            cmd += " -strand_g %s/%s.study.strand" % (self.studyDir, chrName)
        cmd += " -Ne %i" % self.Ne
        cmd += (" -int %i" % startCoord).replace("+","")
        cmd += (" %i" % endCoord).replace("+","")
        cmd += " -o %s" % outF
        cmd += " 2>&1 > %s" % tmpF
        if self.verbose > 1 and self.debug:
            print cmd
            
        exitStatus = os.system( cmd )
        if exitStatus != 0:
            msg = "ERROR: IMPUTE2 returned %i" % exitStatus
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
        if os.path.exists( tmpF ):
            os.remove( tmpF )
        if self.clean:
            for suffix in ["info", "info_by_sample", "summary", "warnings"]:
                fileName = "%s_%s" % (outF, suffix)
                if os.path.exists( fileName ):
                    os.remove( fileName )
        endTime = time.time()
        if self.verbose > 0:
            runLength = datetime.timedelta(seconds=
                                           math.floor(endTime - beginTime))
            msg = " done (%s)\n" % str(runLength)
            sys.stdout.write(msg)
        return outF
    
    
    def launchChrImputation(self, chrName):
        """
        For a given chromosome, define chunk ranges and 
        launch IMPUTE2 on each of them.
        """
        lOutChunkFiles = []
        lCoordRanges = self.getListCoordRangesForEachChunk(chrName, self.chunkLength)
        if self.verbose > 0:
            minCoord = lCoordRanges[0][0]
            maxCoord = lCoordRanges[-1][1]
            msg = "impute '%s'" % chrName
            msg += " (%i->%i,"  % (minCoord, maxCoord)
            msg += " %i chunks)..." % len(lCoordRanges)
            print msg; sys.stdout.flush()
        chunkId = 0
        for chunkCoord in lCoordRanges:
            chunkId += 1
            startCoord = chunkCoord[0]
            endCoord = chunkCoord[1]
            outChunkFile = self.launchChunkImputation(chrName, chunkId,
                                                      startCoord, endCoord)
            lOutChunkFiles.append( outChunkFile )
            if self.debug and chunkId >= 2:
                break
        return lOutChunkFiles
    
    
    def catChunkOutputs(self, chrName, lOutChunkFiles):
        outFile = "%s_chunkAll.impute2" % chrName
        if os.path.exists( outFile ):
            os.remove( outFile )
        if self.verbose > 0:
            print "cat chunks of '%s' into '%s'..." % (chrName, outFile)
            sys.stdout.flush()
        cmd = "cat"
        for outF in lOutChunkFiles:
            if not os.path.exists( outF ):
                msg = "ERROR: chunk output '%s' is missing" % outF
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            cmd += " %s" % outF
        cmd += " > %s" % outFile
        exitStatus = os.system( cmd )
        if exitStatus != 0:
            msg = "ERROR when concatenating the output files for '%s'" % chrName
            sys.stderr.write("%s\n" % msg)
            sys.stderr.write("%s\n" % cmd)
            sys.exit(1)
        if self.gzipOut:
            cmd = "gzip %s" % outFile
            exitStatus = os.system( cmd )
            if exitStatus != 0:
                msg = "ERROR when gzipping the output file for '%s'" % chrName
                sys.stderr.write("%s\n" % msg)
                sys.stderr.write("%s\n" % cmd)
                sys.exit(1)
        if self.clean:
            for outF in lOutChunkFiles:
                os.remove( outF )
                
                
    def run(self):
        self.checkAttributes()
        
        if self.verbose > 0:
            msg = "START %s" % time.strftime("%Y-%m-%d %H:%M:%S")
            startTime = time.time()
            print msg; sys.stdout.flush()
            
        self.getTemplateFiles()
        
        for chrNb in self.lChrNbs:
            chrName = "chr%i" % chrNb
            
            lOutChunkFiles = self.launchChrImputation(chrName)
            
            self.catChunkOutputs(chrName, lOutChunkFiles)
            
        if self.verbose > 0:
            msg = "END %s" % time.strftime("%Y-%m-%d %H:%M:%S")
            endTime = time.time()
            runLength = datetime.timedelta(seconds=
                                           math.floor(endTime - startTime))
            msg += " (%s)" % str(runLength)
            print msg; sys.stdout.flush()
            
            
if __name__ == "__main__":
    i = LaunchImpute2()
    i.setAttributesFromCmdLine()
    i.run()
