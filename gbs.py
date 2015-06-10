#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: perform the bioinformatics aspects of genotyping-by-sequencing
# Copyright (C) 2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

# TODO:
# - check external tools are in PATH
# - add option to ignore R2 files

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
import shlex
import warnings
import shutil
import glob
import random
import string
# import sqlite3
# import numpy as np
# import scipy as sp

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "1.0.1" # http://semver.org/


class TimUtils(object):
    
    def __init__(self):
        pass
    
    @staticmethod
    def user_input(msg):
        if sys.version_info[0] == 2:
            return raw_input(msg)
        elif sys.version_info[0] == 3:
            return input(msg)
        else:
            msg = "Python's major version should be 2 or 3"
            raise ValueError(msg)
        
    @staticmethod
    def uniq_alphanum(length):
        return "".join(random.choice(string.letters+string.digits) \
                       for i in xrange(length))
    
    
class Job(object):
    
    def __init__(self, groupId, cmd, name):
        self.groupId = groupId
        self.cmd = cmd
        self.name = name
        self.queue = None # set by JobGroup upon insertion
        self.duration = None # set by JobGroup upon insertion
        self.memory = None # set by JobGroup upon insertion
        self.id = None # set right after submission
        self.node = None
        
    def submit(self):
        cmd = "echo '%s'" % self.cmd
        cmd += " | qsub -cwd -j y -V -q %s -N %s" % (self.queue, self.name)
        if self.duration:
            pass
        if self.memory:
            pass
        p = Popen(cmd, shell=True, stdout=PIPE).communicate()
        p = p[0].split()[2]
        self.id = int(p)
        return self.id
        
        
class JobGroup(object):
    
    def __init__(self, groupId, scheduler, queue):
        self.id = groupId
        self.scheduler = scheduler
        self.checkScheduler()
        self.queue = queue
        self.checkQueue()
        self.lJobs = [] # filled via self.insert()
        self.lJobIds = [] # filled via self.submit()
        
    def checkScheduler(self):
        if self.scheduler not in ["SGE"]:
            msg = "unknown scheduler '%s'" % self.scheduler
            raise ValueError(msg)
        
    def checkQueue(self):
        if self.scheduler == "SGE":
            p = Popen(["qconf", "-sql"], shell=False, stdout=PIPE).communicate()
            p = p[0].split("\n")
            if self.queue not in p:
                msg = "unknown queue '%s'" % self.queue
                raise ValueError(msg)
            
    def insert(self, iJob):
        self.lJobs.append(iJob)
        self.lJobs[-1].scheduler = self.scheduler
        self.lJobs[-1].queue = self.queue
        
    def submit(self):
        for i in range(len(self.lJobs)):
            jobId = self.lJobs[i].submit()
            self.lJobIds.append(jobId)
            
    def getUnfinishedJobIds(self):
        lUnfinishedJobIds = []
        cmd = "qstat -u '%s'" % os.getlogin()
        cmd += " | grep '%s'" % self.queue
        cmd += " | awk '{print $1}'"
        # print(cmd) # debug
        p = Popen(cmd, shell=True, stdout=PIPE).communicate()
        # p1 = Popen(shlex.split("qstat -u '%s'" % os.getlogin()), stdout=PIPE)
        # p2 = Popen(shlex.split("grep '%s'" % self.queue), stdin=p1.stdout,
        #            stdout=PIPE)
        # p3 = Popen(shlex.split("awk '{print $1}'"), stdin=p2.stdout,
        #            stdout=PIPE)
        # p = p3.communicate()
        # print(p) # debug
        p = p[0].split("\n")[:-1]
        lUnfinishedJobIds = [int(jobId) for jobId in p
                             if int(jobId) in self.lJobIds]
        return lUnfinishedJobIds
    
    def wait(self):
        time.sleep(2)
        if len(self.getUnfinishedJobIds()) == 0:
            return
        time.sleep(5)
        if len(self.getUnfinishedJobIds()) == 0:
            return
        time.sleep(10)
        if len(self.getUnfinishedJobIds()) == 0:
            return
        time.sleep(30)
        if len(self.getUnfinishedJobIds()) == 0:
            return
        while True:
            time.sleep(60)
            if len(self.getUnfinishedJobIds()) == 0:
                return
            
            
class JobManager(object):
    
    def __init__(self, projectId):
        self.projectId = projectId
        self.dbPath = ""
        self.groupId2group = {}
        
    def __getitem__(self, jobGroupId):
        return self.groupId2group[jobGroupId]
    
    def insert(self, iJobGroup):
        self.groupId2group[iJobGroup.id] = iJobGroup
        
        
class GbsSample(object):
    """
    Corresponds to some DNA present on a given lane. This DNA was extracted 
    from an individual which can also have some DNA on other lanes, which will 
    then be considered as other "samples".
    """
    
    def __init__(self, id):
        self.id = id
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
        self.dDemultiplexedFastqFiles = {}
        
    def __str__(self):
        txt = "id=%s" % self.id
        txt += ";individual=%s" % self.individual
        txt += ";generation=%s" % self.generation
        txt += ";species=%s" % self.species
        txt += ";barcode=%s" % self.barcode
        txt += ";seqCenter=%s" % self.seqCenter
        txt += ";seqPlatform=%s" % self.seqPlatform
        txt += ";seqPlatformModel=%s" % self.seqPlatformModel
        txt += ";flowcell=%s" % self.flowcell
        txt += ";lane=%s" % self.lane
        txt += ";date=%s" % self.date
        txt += ";initFastqFile1=%s" %self.initFastqFile1
        txt += ";initFastqFile2=%s" %self.initFastqFile2
        return txt
    
    def setDemultiplexedFastqFiles(self, pathToDir):
        lFilesR1 = glob.glob("%s/*_%s*_R1.fastq.gz" % (pathToDir,
                                                       self.individual))
        if len(lFilesR1) != 1:
            msg = "can't find R1 file for sample '%s' in '%s'" % (self.id,
                                                                  pathToDir)
            raise ValueError(msg)
        self.dDemultiplexedFastqFiles["R1"] = lFilesR1[0]
        lFilesR2 = glob.glob("%s/*_%s*_R2.fastq.gz" % (pathToDir,
                                                       self.individual))
        if len(lFilesR2) == 1:
            self.dDemultiplexedFastqFiles["R2"] = lFilesR2[0]
            
    def clean(self, adpFwd, adpRev, outDir, iJobGroup):
        cmd = "cutadapt"
        cmd += " -a %s" % adpFwd
        cmd += " -A %s" % adpRev
        cmd += " -o %s_clean_R1.fastq.gz" % self.individual
        cmd += " -p %s_clean_R2.fastq.gz" % self.individual
        cmd += " -e 0.2"
        cmd += " -O 3"
        cmd += " -m 35"
        cmd += " -U 3"
        cmd += " -q 20,20"
        cmd += " --max-n 0.2"
        # cmd += " --maximum-length 150"
        cmd += " %s" % self.dDemultiplexedFastqFiles["R1"]
        if "R2" in self.dDemultiplexedFastqFiles:
            cmd += " %s" % self.dDemultiplexedFastqFiles["R2"]
        print(cmd)
        sys.exit(1)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = Job(iJobGroup.id, cmd, jobName)
        iJobGroup.insert(iJob)
        
        
class GbsLane(object):
    
    def __init__(self, laneId, flowcell, number):
        self.id = laneId # flowcell identifier + "_" + lane number
        self.flowcell = flowcell
        self.number = number
        self.dSamples = {}
        self.dInitFastqFileSymlinks = {} # key(s): R1 (and R2, optional)
        
    def insert(self, iSample):
        if iSample.id in self.dSamples:
            msg = "sample '%s' already in lane '%s'" % (iSample.id, self.id)
            raise ValueError(msg)
        self.dSamples[iSample.id] = iSample
        
    def getInitFastqFiles(self):
        lFiles = []
        for sampleId,iSample in self.dSamples.items():
            if iSample.initFastqFile1 not in lFiles:
                lFiles.append(iSample.initFastqFile1)
            if iSample.initFastqFile2 and iSample.initFastqFile2 not in lFiles:
                lFiles.append(iSample.initFastqFile2)
        if len(lFiles) == 0:
            msg = "lane '%s' has no initial fastq file" % self.id
            raise ValueError(msg)
        if len(lFiles) > 2:
            msg = "lane '%s' has more than 2 initial fastq files" % self.id
            raise ValueError(msg)
        return lFiles
    
    def setInitFastqFileSymlinks(self, lInitFastqFiles, pathToDir):
        self.dInitFastqFileSymlinks["R1"] = [lInitFastqFiles[0],
                                             "%s/%s_R1.fastq.gz" \
                                             % (pathToDir, self.id)]
        if len(lInitFastqFiles) > 1:
            self.dInitFastqFileSymlinks["R2"] = [lInitFastqFiles[1],
                                                 "%s/%s_R2.fastq.gz" \
                                                 % (pathToDir, self.id)]
            
    def initQc(self, outDir, iJobGroup):
        for Ri,lFiles in self.dInitFastqFileSymlinks.items():
            cmd = "time fastqc -o %s %s" % (outDir, lFiles[1])
            jobName = "stdout_%s_%s_%s" % (iJobGroup.id, self.id, Ri)
            iJob = Job(iJobGroup.id, cmd, jobName)
            iJobGroup.insert(iJob)
            
    def saveBarcodeFile(self, format="fasta"):
        fileName = "barcodes_%s.fa" % self.id
        fileHandle = open(fileName, "w")
        for sample,iSample in self.dSamples.items():
            fileHandle.write(">%s\n%s\n" % (iSample.individual,
                                            iSample.barcode))
        fileHandle.close()
        return fileName
        
    def demultiplex(self, outDir, enzyme, iJobGroup):
        if len(self.dInitFastqFileSymlinks) < 2:
            msg = "can't demultiplex (yet) lane '%s' if single-end" % self.id
            raise ValueError(msg)
        cmd = "demultiplex.py"
        cmd += " --idir %s" % os.path.dirname(self.dInitFastqFileSymlinks["R1"][1])
        cmd += " --ifq1 %s" % os.path.basename(self.dInitFastqFileSymlinks["R1"][1])
        cmd += " --ifq2 %s" % os.path.basename(self.dInitFastqFileSymlinks["R2"][1])
        cmd += " --it %s" % self.saveBarcodeFile(self)
        cmd += " --ofqp %s/%s" % (outDir, self.id)
        cmd += " --met %s" % "4c"
        cmd += " --re %s" % enzyme
        cmd += " --chim 1"
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = Job(iJobGroup.id, cmd, jobName)
        iJobGroup.insert(iJob)
        
    def setDemultiplexedFastqFiles(self, pathToDir):
        for sampleId,iSample in self.dSamples.items():
            iSample.setDemultiplexedFastqFiles("%s/%s" % (pathToDir, self.id))
            
    def clean(self, adpFwd, adpRev, outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.clean(adpFwd, adpRev, outDir, iJobGroup)
            
            
class Gbs(object):
    
    def __init__(self):
        self.verbose = 1
        self.lSteps = []
        self.infoFile = ""
        self.dictFile = ""
        self.scheduler = "SGE"
        self.queue = "normal.q"
        self.projectId = ""
        self.enzyme = "ApeKI"
        self.jobManager = None
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
        self.dLanes = {}
        self.ind2samples = {}
        self.lDirSteps = ["lane_qualities",
                          "demultiplexing_lanes",
                          "cleaning_samples",
                          "mapping_samples",
                          "realign_bqsr",
                          "calling_variants"]
        self.adpFile = ""
        self.adapters = {}
        
        
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
        msg += "      --proj\tname of the project\n"
        msg += "      --schdlr\tname of the cluster scheduler (default=SGE)\n"
        msg += "      --queue\tname of the cluster queue (default=normal.q)\n"
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
        msg += "      --adp\tpath to the file containing the adapters\n"
        msg += "\t\tsame format as FastQC: name<tab>sequence\n"
        msg += "\t\tname: at least 'adpFwd' (also 'adpRev' if paired-end)\n"
        msg += "\t\tsequence: from 5' to 3'\n"
        msg += "      --dict\tpath to the 'dict' file (SAM header with @SQ tags)\n"
        msg += "      --enz\tname of the restriction enzyme (default=ApeKI)\n"
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
                                         "proj=", "step=", "info=", "dict=",
                                         "schdlr=", "queue=", "enz=", "adp="])
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
            elif o == "--proj":
                self.projectId = a
            elif o == "--schdlr":
                self.scheduler = a
            elif o == "--queue":
                self.queue = a
            elif o == "--step":
                self.lSteps = a.split("-")
            elif o == "--info":
                 self.infoFile = a
            elif o == "--dict":
                 self.dictFile = a
            elif o == "--enz":
                self.enzyme = a
            elif o == "--adp":
                self.adpFile = a
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if self.projectId == "":
            msg = "ERROR: missing compulsory option --proj"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "_" in self.projectId:
            msg = "ERROR: forbidden underscore '_' in project identifier '%s'" \
                  % self.projectId
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
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
        if self.scheduler == "":
            msg = "ERROR: missing compulsory option --schdlr"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.scheduler == "OGE":
            self.scheduler = "SGE"
        if self.queue == "":
            msg = "ERROR: missing compulsory option --queue"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        self.jobManager = JobManager(self.projectId)
        if self.lSteps == []:
            msg = "ERROR: missing compulsory option --step"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        for i in range(len(self.lDirSteps)):
            self.lDirSteps[i] = "%s/%s_%s" % (os.getcwd(), self.projectId,
                                              self.lDirSteps[i])
        if "3" in self.lSteps:
            if self.adpFile == "":
                msg = "ERROR: missing compulsory option --adp"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not os.path.exists(self.adpFile):
                msg = "ERROR: can't find file %s" % self.adpFile
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
                
                
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
            
            # create and fill a "GbsSample" object
            individual = tokens[self.infoCol2idx["individual"]]
            flowcell = tokens[self.infoCol2idx["flowcell"]]
            laneNum = int(tokens[self.infoCol2idx["lane"]])
            sampleId = "%s_%s_%i" % (individual, flowcell, laneNum)
            iSample = GbsSample(sampleId)
            iSample.individual = individual
            iSample.flowcell = flowcell
            iSample.lane = laneNum
            iSample.generation = tokens[self.infoCol2idx["generation"]]
            iSample.species = tokens[self.infoCol2idx["species"]]
            iSample.barcode = tokens[self.infoCol2idx["barcode"]]
            iSample.seqCenter = tokens[self.infoCol2idx["seq_center"]]
            iSample.seqPlatform = tokens[self.infoCol2idx["seq_platform"]]
            iSample.seqPlatformModel = tokens[self.infoCol2idx["seq_platform_model"]]
            iSample.date = tokens[self.infoCol2idx["date"]]
            if not os.path.exists(tokens[self.infoCol2idx["fastq_file_R1"]]):
                msg = "file '%s' doesn't exist" \
                      % tokens[self.infoCol2idx["fastq_file_R1"]]
                raise OSError(msg)
            if not tokens[self.infoCol2idx["fastq_file_R1"]].endswith(".gz"):
                msg = "file '%s' should be gzipped" \
                      % tokens[self.infoCol2idx["fastq_file_R1"]]
                raise ValueError(msg)
            iSample.initFastqFile1 = tokens[self.infoCol2idx["fastq_file_R1"]]
            if tokens[self.infoCol2idx["fastq_file_R2"]] != "":
                if not os.path.exists(tokens[self.infoCol2idx["fastq_file_R2"]]):
                    msg = "file '%s' doesn't exist" \
                          % tokens[self.infoCol2idx["fastq_file_R2"]]
                    raise OSError(msg)
                if not tokens[self.infoCol2idx["fastq_file_R2"]].endswith(".gz"):
                    msg = "file '%s' should be gzipped" \
                          % tokens[self.infoCol2idx["fastq_file_R2"]]
                    raise ValueError(msg)
                iSample.initFastqFile2 = tokens[self.infoCol2idx["fastq_file_R2"]]
            self.dSamples[iSample.id] = iSample
            
            laneId = "%s_%i" % (flowcell, laneNum)
            if laneId not in self.dLanes:
                self.dLanes[laneId] = GbsLane(laneId, flowcell, laneNum)
            self.dLanes[laneId].insert(iSample)
            
            if individual not in self.ind2samples:
                self.ind2samples[individual] = []
            self.ind2samples[individual].append(sampleId)
            
            
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
            msg += "\nnb of lanes: %i" % len(self.dLanes)
            msg += "\nnb of individuals: %i" % len(self.ind2samples)
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()


    def enterStepDir(self, dirStep):
        if os.path.isdir(dirStep):
            msg = "directory '%s' already exists" % dirStep
            warnings.warn(msg, UserWarning)
            wantRmvDir = TimUtils.user_input("Do you want to remove the directory '%s'? [y/n] " % dirStep)
            if wantRmvDir == "y":
                shutil.rmtree(dirStep)
            else:
                raise OSError("can't continue")
        os.mkdir(dirStep)
        cwd = os.getcwd()
        os.chdir(dirStep)
        return cwd


    def beginStep(self, stepNum):
        if self.verbose > 0:
            msg = "perform step %i ..." % stepNum
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        cwd = self.enterStepDir(self.lDirSteps[stepNum - 1])
        return cwd
    
    
    def setPathsToInputFiles(self, stepNum):
        if stepNum in [1, 2]:
            for laneId,iLane in self.dLanes.items():
                lFiles = iLane.getInitFastqFiles()
                iLane.setInitFastqFileSymlinks(lFiles, self.lDirSteps[0])
        if stepNum == 3:
            for laneId,iLane in self.dLanes.items():
                iLane.setDemultiplexedFastqFiles(self.lDirSteps[1])
                
                
    def endStep(self, stepNum, cwd):
        os.chdir(cwd)
        if self.verbose > 0:
            msg = "step %i is done" % stepNum
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def makeSymlinksToInputFastqFiles(self):
        """
        Make one symlink per lane for R1 (and R2 if necessary).
        """
        for lane,iLane in self.dLanes.items():
            for r,lFiles in iLane.dInitFastqFileSymlinks.items():
                os.symlink(lFiles[0], lFiles[1])
                
                
    def makeGroupJobId(self, prefix):
        return "%s_%s" % (prefix, TimUtils.uniq_alphanum(5))
    
    
    def launchFastqcOnInputFastqFiles(self):
        if self.verbose > 0:
            msg = "assess quality per lane ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step1-fastqc" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            iLane.initQc(self.lDirSteps[0], iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all quality jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def combineFastqcResults(self):
        # if self.verbose > 0:
        #     msg = "combine FastQC results ..."
        #     sys.stdout.write("%s\n" % msg)
        #     sys.stdout.flush()
        pass
    
    
    def step1(self):
        self.setPathsToInputFiles(1)
        cwd = self.beginStep(1)
        self.makeSymlinksToInputFastqFiles()
        self.launchFastqcOnInputFastqFiles()
        self.combineFastqcResults()
        self.endStep(1, cwd)
        
        
    def loadAdpFile(self):
        if self.verbose > 0:
            msg = "load adapter file ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        adpHandle = open(self.adpFile)
        while True:
            line = adpHandle.readline()
            if line == "":
                break
            tokens = line.rstrip("\n").split("\t")
            if len(tokens) != 2:
                msg = "adapter file '%s' should have two columns separated by a tabulation" \
                      % self.adpFile
                raise ValueError(msg)
            self.adapters[tokens[0]] = tokens[1]
        adpHandle.close()
        if "adpFwd" not in self.adapters or "adpRev" not in self.adapters:
            msg = "'adpFwd' or 'adpRev' missing from adapter file '%s'" % \
                  self.adpFile
            raise ValueError(msg)
        
        if self.verbose > 0:
            msg = "nb of adapters: %i" % len(self.adapters)
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def demultiplexInputFastqFiles(self):
        if self.verbose > 0:
            msg = "demultiplex samples per lane ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step2-demultiplex" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[1], laneId)
            os.mkdir(outDir)
            iLane.demultiplex(outDir, self.enzyme, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all demultiplexing jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step2(self):
        self.setPathsToInputFiles(2)
        cwd = self.beginStep(2)
        self.demultiplexInputFastqFiles()
        self.endStep(2, cwd)
        
        
    def cleanDemultiplexedFiles(self):
        if self.verbose > 0:
            msg = "clean reads per sample ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step3-clean" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[2], laneId)
            os.mkdir(outDir)
            iLane.clean(self.adapters["adpFwd"], self.adapters["adpRev"],
                        outDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all cleaning jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step3(self):
        self.loadAdpFile()
        self.setPathsToInputFiles(3)
        cwd = self.beginStep(3)
        self.cleanDemultiplexedFiles()
        self.endStep(3, cwd)
        
        
    def mapCleanedReads(self):
        if self.verbose > 0:
            msg = "map reads per sample ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step4-map" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for sampleId,iSample in self.dSamples.items():
            iSample.map(iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all mapping jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step4(self):
        cwd = self.beginStep(4)
        self.mapCleanedReads()
        self.endStep(4, cwd)
        
        
    def localRealignment(self):
        pass
    
    def baseQualityScoreRecalibration(self):
        pass
    
    
    def step5(self):
        cwd = self.beginStep(5)
        self.localRealignment()
        self.baseQualityScoreRecalibration()
        self.endStep(5, cwd)
        
        
    def variantCallingPerIndividual(self):
        pass
    
    def jointGenotyping(self):
        pass

    def variantRecalibration(self):
        pass
    
    
    def step6(self):
        cwd = self.beginStep(6)
        self.variantCallingPerIndividual()
        self.jointGenotyping()
        self.variantRecalibration()
        self.endStep(6, cwd)
        
        
    def run(self):
        self.loadInfoFile()
        
        if "1" in self.lSteps:
            # QC: parallelize over lanes
            self.step1()
            
        if "2" in self.lSteps:
            # demultiplex: parallelize over lanes
            self.step2()
            
        if "3" in self.lSteps:
            # clean: parallelize over samples
            self.step3()
            
        if "4" in self.lSteps:
            # map: parallelize over samples
            self.step4()
            
        if "5" in self.lSteps:
            # recalibrate: parallelize over lanes
            self.step5()
            
        if "6" in self.lSteps:
            # genotype: parallelize over individuals
            self.step6()
            
            
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
