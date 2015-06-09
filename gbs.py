#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: perform the bioinformatics aspects of genotyping-by-sequencing
# Copyright (C) 2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import sys
import os
import getopt
import time
import datetime
from subprocess import Popen, PIPE
import math
import gzip
# import numpy as np
# import scipy as sp

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "1.0.0" # http://semver.org/


class GbsSample(object):
    def __init__(self, name):
        self.name = name
        self.individual = ""
        self.generation = ""
        self.species = ""
        self.barcode = ""
        self.seqCenter = ""
        self.seqPlatform = ""
        self.seqPlatformModel = ""
        self.flowcell = ""
        self.lane = ""
        self.date = ""
        self.initFastqFile1 = ""
        self.initFastqFile2 = None
        
        
class Gbs(object):
    
    def __init__(self):
        self.verbose = 1
        self.lSteps = []
        self.infoFile = ""
        self.dictFile = ""
        self.clusterQueue = ""
        self.projectDir = ""
        self.infoCol2idx = {"individual": None,
                            "generation": None,
                            "species": None,
                            "barcode": None,
                            "seq_center": None,
                            "seq_platform": None,
                            "seq_platform_model": None,
                            "flowcell": None,
                            "lane": None,
                            "date": None,
                            "fastq_file_R1": None,
                            "fastq_file_R2": None}
        self.dSamples = {}
        self.laneId2samples = {}
        self.ind2samples = {}
        self.dirStep1 = "lane_qualities"
        self.dirStep2 = "demultiplexing_lanes"
        self.dirStep3 = "cleaning_samples"
        self.dirStep4 = "mapping_samples"
        self.dirStep5 = "realign_bqsr"
        self.dirStep6 = "calling_variants"
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' performs the bioinformatics aspects of genotyping-by-sequencing.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "      --step\tstep(s) to perform (1/2/3/4, can be 1-2)\n"
        msg += "\t\t1: raw read quality (with FastQC v >= 0.11.2)\n"
        msg += "\t\t2: demultiplexing\n"
        msg += "\t\t3: cleaning (with CutAdapt v >= 1.8)\n"
        msg += "\t\t4: mapping (with BWA MEM v >= 0.7.12 and Samtools v >= 1.1)\n"
        msg += "\t\t5: local realignment and BQSR (with GATK v >= 3.3)\n"
        msg += "\t\t6: SNP and genotype calling (with GATK HC and/or BcfTools)\n"
        msg += "      --info\tpath to the 'information' file\n"
        msg += "\t\tthe file should be encoded in ASCII\n"
        msg += "\t\tthe first row should be a header with column names\n"
        msg += "\t\teach sample should have one and only one row\n"
        msg += "\t\tany two columns should be separated with one tabulation\n"
        msg += "\t\tcolumns can be in any order\n"
        msg += "\t\t12 columns are compulsory:\n"
        msg += "\t\t individual (e.g. 'Col-0')\n"
        msg += "\t\t generation (empty or 0/1/etc, with 0 for founders)\n"
        msg += "\t\t species (e.g. 'Arabidopsis thaliana')\n"
        msg += "\t\t barcode (e.g. 'ATGG')\n"
        msg += "\t\t seq_center (e.g. 'Broad Institute', 'GenoToul', etc)\n"
        msg += "\t\t seq_platform (e.g. 'ILLUMINA', see SAM format specification)\n"
        msg += "\t\t seq_platform_model (e.g. 'HiSeq 2000')\n"
        msg += "\t\t flowcell (e.g. 'C5YMDACXX')\n"
        msg += "\t\t lane (e.g. '3')\n"
        msg += "\t\t date (e.g. '2015-01-15', see SAM format specification)\n"
        msg += "\t\t fastq_file_R1 (one per lane, gzip-compressed)\n"
        msg += "\t\t fastq_file_R2 (one per lane, gzip-compressed)\n"
        msg += "      --dict\tpath to the 'dict' file (SAM header with @SQ tags)\n"
        msg += "      --queue\tname of the cluster queue\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s --step 1 --info info_gbs.txt\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@supagro.inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        
        The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012).
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2015 Institut National de la Recherche Agronomique.\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:i:",
                                        ["help", "version", "verbose=",
                                         "step=", "info=", "dict="])
        except getopt.GetoptError as err:
            sys.stderr.write("%s\n\n" % str(err))
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
            elif o == "--step":
                self.lSteps = a.split("-")
            elif o == "--info":
                 self.infoFile = a
            elif o == "--dict":
                 self.dictFile = a
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if self.lSteps == []:
            msg = "ERROR: missing compulsory option --step"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "1" in self.lSteps:
            if self.infoFile == "":
                msg = "ERROR: missing compulsory option --info"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not os.path.exists(self.infoFile):
                msg = "ERROR: can't find file %s" % self.infoFile
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
        self.projectDir = os.getcwd()
        
        
    def loadHeaderInfoFile(self, line):
        """
        Set values to self.infoCol2idx.
        """
        try:
            line.decode('ascii')
        except UnicodeDecodeError as err:
            raise
        tokens = line.rstrip("\n").split("\t")
        if len(tokens) < 12:
            msg = "header should have at least 12 tab-separated columns"
            raise ValueError(msg)
        for idx,tok in enumerate(tokens):
            if tok in self.infoCol2idx:
                self.infoCol2idx[tok] = idx
        for infoCol,idx in self.infoCol2idx.items():
            if idx is None:
                msg = "column '%s' not found in information file" % infoCol
                raise ValueError(msg)
            
            
    def loadContentInfoFile(self, lines):
        for line in lines:
            tokens = line.rstrip("\n").split("\t")
            
            individual = tokens[self.infoCol2idx["individual"]]
            flowcell = tokens[self.infoCol2idx["flowcell"]]
            lane = tokens[self.infoCol2idx["lane"]]
            name = "%s_%s_%s" % (individual, flowcell, lane)
            self.dSamples[name] = GbsSample(name)
            self.dSamples[name].individual = individual
            self.dSamples[name].flowcell = flowcell
            self.dSamples[name].lane = lane
            self.dSamples[name].name = name
            laneId = "%s_%s" % (flowcell, lane)
            if laneId not in self.laneId2samples:
                self.laneId2samples[laneId] = []
            self.laneId2samples[laneId].append(name)
            if individual not in self.ind2samples:
                self.ind2samples[individual] = []
            self.ind2samples[individual].append(name)
            
            self.dSamples[name].generation \
                = tokens[self.infoCol2idx["generation"]]
            self.dSamples[name].species \
                = tokens[self.infoCol2idx["species"]]
            self.dSamples[name].barcode \
                = tokens[self.infoCol2idx["barcode"]]
            self.dSamples[name].seq_center \
                = tokens[self.infoCol2idx["seq_center"]]
            self.dSamples[name].seq_platform \
                = tokens[self.infoCol2idx["seq_platform"]]
            self.dSamples[name].seq_platform_model \
                = tokens[self.infoCol2idx["seq_platform_model"]]
            self.dSamples[name].date \
                = tokens[self.infoCol2idx["date"]]
            if not os.path.exists(tokens[self.infoCol2idx["fastq_file_R1"]]):
                msg = "file '%s' doesn't exist" \
                      % tokens[self.infoCol2idx["fastq_file_R1"]]
                raise OSError(msg)
            if not tokens[self.infoCol2idx["fastq_file_R1"]].endswith(".gz"):
                msg = "file '%s' should be gzipped" \
                      % tokens[self.infoCol2idx["fastq_file_R1"]]
                raise ValueError(msg)
            self.dSamples[name].fastq_file_R1 \
                = tokens[self.infoCol2idx["fastq_file_R1"]]
            if tokens[self.infoCol2idx["fastq_file_R2"]] != "":
                if not os.path.exists(tokens[self.infoCol2idx["fastq_file_R2"]]):
                    msg = "file '%s' doesn't exist" \
                          % tokens[self.infoCol2idx["fastq_file_R2"]]
                    raise OSError(msg)
                if not tokens[self.infoCol2idx["fastq_file_R2"]].endswith(".gz"):
                    msg = "file '%s' should be gzipped" \
                          % tokens[self.infoCol2idx["fastq_file_R2"]]
                    raise ValueError(msg)
                self.dSamples[name].fastq_file_R2 \
                    = tokens[self.infoCol2idx["fastq_file_R2"]]
                
                
    def loadInfoFile(self):
        if self.verbose > 0:
            msg = "load information file '%s' ..." % self.infoFile
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        infoHandle = open(self.infoFile)

        lines = infoHandle.readlines()
        self.loadHeaderInfoFile(lines[0])
        self.loadContentInfoFile(lines[1:])
        
        infoHandle.close()
        if self.verbose > 0:
            msg = "nb of samples: %i" % len(self.dSamples)
            msg += "\nnb of lanes: %i" % len(self.laneId2samples)
            msg += "\nnb of individuals: %i" % len(self.ind2samples)
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()


    def enterStepDir(self, dirStep):
        if os.path.isdir(dirStep):
            msg = "directory '%s' already exists" % dirStep
            raise OSError(msg)
        os.mkdir(dirStep)
        os.chdir(dirStep)
        
        
    def makeSymlinksToInputFastqFiles(self):
        """
        Make one symlink per lane for R1 (and R2 if necessary).
        """
        for lane,lSamples in self.laneId2samples.items():
            iSample = self.dSamples[lSamples[0]]
            bn1 = os.path.basename(iSample.fastq_file_R1)
            os.symlink(iSample.fastq_file_R1, "%s/%s/%s_R1.fastq.gz" \
                       % (self.projectDir, self.dirStep1, lane))
            if iSample.fastq_file_R2 is not None:
                bn2 = os.path.basename(iSample.fastq_file_R2)
                os.symlink(iSample.fastq_file_R2, "%s/%s/%s_R2.fastq.gz" \
                           % (self.projectDir, self.dirStep1, lane))
                
                
    def makeGroupJobId(self):
        pass


    def submitJob(self, cmd, groupJobId):
        pass
        
        
    def launchFastqcOnInputFastqFiles(self):
        if self.verbose > 0:
            msg = "launch FastQC on original fastq files ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        groupJobId = self.makeGroupJobId()
        for lane in self.laneId2samples:
            for sample in self.laneId2samples[lane]:
                iSample = self.dSamples[sample]
                cmd = "time fastqc -o %s/%s %s" % \
                      (self.projectDir, self.dirStep1,
                       iSample.initFastqFile1)
                self.submitJob(cmd, groupJobId)
                if iSample.initFastqFile2:
                    cmd = "time fastqc -o %s/%s %s" % \
                          (self.projectDir, self.dirStep1,
                           iSample.initFastqFile1)
                    self.submitJob(cmd, groupJobId)
                    
                    
    def combineFastqcResults(self):
        if self.verbose > 0:
            msg = "combine FastQC results ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        pass
    
    
    def step1(self):
        step = "1"
        if self.verbose > 0:
            msg = "perform step %s ..." % step
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        self.enterStepDir(self.dirStep1)
        self.makeSymlinksToInputFastqFiles()
        self.launchFastqcOnInputFastqFiles()
        self.combineFastqcResults()
        os.chdir(self.projectDir)
        
        if self.verbose > 0:
            msg = "step %s is done" % step
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def demultiplexInputFastqFiles(self):
        pass
    
    
    def step2(self):
        step = "2"
        if self.verbose > 0:
            msg = "perform step %s ..." % step
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        self.enterStepDir(self.dirStep2)
        self.demultiplexInputFastqFiles()
        os.chdir(self.projectDir)
        
        if self.verbose > 0:
            msg = "step %s is done" % step
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def cleanDemultiplexedFiles(self):
        pass
    
    
    def step3(self):
        step = "3"
        if self.verbose > 0:
            msg = "perform step %s ..." % step
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        self.enterStepDir(self.dirStep3)
        self.cleanDemultiplexedFiles()
        os.chdir(self.projectDir)
        
        if self.verbose > 0:
            msg = "step %s is done" % step
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def mapCleanedReads(self):
        pass
    
    
    def step4(self):
        step = "4"
        if self.verbose > 0:
            msg = "perform step %s ..." % step
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        self.enterStepDir(self.dirStep4)
        self.mapCleanedReads()
        os.chdir(self.projectDir)
        
        if self.verbose > 0:
            msg = "step %s is done" % step
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
    
    
    def localRealignment(self):
        pass
    
    def baseQualityScoreRecalibration(self):
        pass
    
    def variantCallingPerIndividual(self):
        pass
    
    def jointGenotyping(self):
        pass

    def variantRecalibration(self):
        pass
    
    
    def run(self):
        if "1" in self.lSteps:
            self.loadInfoFile()
            self.step1()
            
        if "2" in self.lSteps:
            self.step2()
            
        if "3" in self.lSteps:
            self.step3()
            
        if "4" in self.lSteps:
            self.step4()
            
        if "5" in self.lSteps:
            # create dir "realign_bqsr"
            # launch GATK
            self.localRealignment()
            self.baseQualityScoreRecalibration()
            
        if "6" in self.lSteps:
            # create dir "calling_variants"
            # launch GATK HaplotypeCaller
            self.variantCallingPerIndividual()
            self.jointGenotyping()
            self.variantRecalibration()
            
            
if __name__ == "__main__":
    i = Gbs()
    
    i.setAttributesFromCmdLine()
    
    i.checkAttributes()
    
    if i.verbose > 0:
        startTime = time.time()
        msg = "START %s %s %s" % (os.path.basename(sys.argv[0]),
                                  progVersion,
                                  time.strftime("%Y-%m-%d %H:%M:%S"))
        msg += "\ncmd-line: %s" % ' '.join(sys.argv)
        msg += "\ncwd: %s" % os.getcwd()
        print(msg); sys.stdout.flush()
        
    i.run()
    
    if i.verbose > 0:
        msg = "END %s %s %s" % (os.path.basename(sys.argv[0]),
                                progVersion,
                                time.strftime("%Y-%m-%d %H:%M:%S"))
        endTime = time.time()
        runLength = datetime.timedelta(seconds=
                                       math.floor(endTime - startTime))
        msg += " (%s" % str(runLength)
        if "linux" in sys.platform:
            p = Popen(["grep", "VmHWM", "/proc/%s/status" % os.getpid()],
                      shell=False, stdout=PIPE).communicate()
            maxMem = p[0].split()[1]
            msg += "; %s kB)" % maxMem
        else:
            msg += ")"
        print(msg); sys.stdout.flush()
